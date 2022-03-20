#include <iostream>
#include <map>
#include <deque>
#include <utility>
#include <vector>

#include "tfm_index.hpp"
extern "C" {
#include "utils.h"
}

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;


struct Dict {
    uint8_t *d;  // pointer to the dictionary
    long *end;   // end[i] is the index of the ending symbol of the i-th phrase
    long dsize;  // dicionary size in symbols
    int dwords;  // the number of phrases of the dicionary
};

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << "w BASENAME" << endl;
	cerr << "BASENAME:" << endl;
	cerr << "  Name prefix of the data files to process." << endl;
};

Dict read_dictionary(char *filename) {
        FILE *g = open_aux_file(filename, EXTDICT,"rb");
        fseek(g, 0, SEEK_END);
        long dsize = ftell(g);
        if(dsize < 0) die("ftell dictionary");
        if(dsize <= 1+4) die("invalid dictionary file");
        cout  << "Dictionary file size: " << dsize << endl;
        #if !M64
        if(dsize > 0x7FFFFFFE) {
            printf("Dictionary size greater than  2^31-2!\n");
            printf("Please use 64 bit version\n");
            exit(1);
        }
        #endif

        uint8_t *d = new uint8_t[dsize];
        rewind(g);
        long e = fread(d, 1, dsize, g);
        if(e!=dsize) die("fread");
        fclose(g);

        int dwords = 0;
        for (int i = 0; i < dsize; i++) {
            if (d[i] == EndOfWord) dwords++;
        }
        cout << "Dictionary contains " << dwords << " words" << endl;

        long *end= new long[dwords];
        int cnt = 0;
        for (int i = 0; i < dsize; i++) {
            if (d[i] == EndOfWord) end[cnt++] = i;
        }

        Dict res = {d, end, dsize, dwords};
        return res;
}

void write_occs(char *basename, Dict dict, vector<vector<uint32_t>> &sources) {
    char *name;
    asprintf(&name, "%s.%s", basename, "occ");
    FILE *focc = fopen(name, "wb");
    for (int i = 0; i < dict.dwords; i++) {
        uint32_t cnt = sources[i].size();
        size_t s = fwrite(&cnt, sizeof(uint32_t), 1, focc);
        if (s != 1) {
            cout << "Error writing to OCC file" << endl;
            exit(1);
        }
    }
    cout << "OCC file successfully written!" << endl;
    fclose(focc);
}

void write_ilist(char *basename, Dict dict, vector<vector<uint32_t>> &sources) {
    char *name;
    asprintf(&name, "%s.%s", basename, "ilist");
    FILE *filist = fopen(name, "wb");
    for (int i = 0; i < dict.dwords; i++) {
        for (int j = 0; j < (int)sources[i].size(); j++) {
            size_t s = fwrite(&sources[i][j], sizeof(uint32_t), 1, filist);
            if (s != 1) {
                cout << "Error writing to ILIST file" << endl;
                exit(1);
            }
        }
    }
    cout << "ILIST file successfully written!" << endl;
    fclose(filist);
}

int main(int argc, char **argv) {
	//check parameters
	if (argc < 2) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}
        int w = atoi(argv[1]);
        // load the dictionary and the tunneled WG
        struct Dict dict = read_dictionary(argv[2]);
        char *name;
        asprintf(&name, "%s.%s", argv[2], "tunnel");
        tfm_index<> tfm;
        load_from_file(tfm, name);
        cout << "TFM loaded" << endl;

        vector<vector<uint32_t>> phrase_sources(dict.dwords+1); // the list of ourgoing edges for each vertex
        vector<vector<uint32_t>> phrase_destinations(dict.dwords+1); // the list of ourgoing edges for each vertex

        for (int i = 0; i < tfm.L.size(); i++) {
            uint32_t act_char = tfm.L[i];
            phrase_sources[act_char].push_back(i);
        }
        cout << "Inverted lists created!" << endl;

        // write the numbers of occurrences for each phrase
        write_occs(argv[2], dict, phrase_sources);
        // write the concatenation of the 'inverted lists' of the phrases
        write_ilist(argv[2], dict, phrase_sources);

        vector<vector<bool>> phrase_din(dict.dwords+1);
        vector<vector<bool>> phrase_dout(dict.dwords+1);
        vector<vector<uint8_t>> phrase_outlabels(dict.dwords+1); // the list of ourgoing edges for each vertex
        uint64_t last_size = 0, bitvector_size = 0;
        for (int c = 0; c < dict.dwords+1; c++) {
            vector<bool> seen(tfm.dout_rank(tfm.L.size()), false); // mark processed vertices
            for (int i = 0; i < (int)phrase_sources[c].size(); i++) {

                uint32_t pos = phrase_sources[c][i];
                auto is = tfm.L.inverse_select(pos).first; // inverse_select returns rank(L[i], i), L[i])
                uint32_t new_i = tfm.C[c] + is; // we want to include the is-th position into the rank
                uint32_t din_rank = tfm.din_rank(new_i + 1);
                if (!seen[din_rank]) {
                    phrase_din[c].push_back(1);
                    phrase_dout[c].push_back(1);
                    seen[din_rank] = true;

                    // process all the outgoing edges
                    int act_i = new_i;
                    while (1) {
                        uint32_t act_phrase = tfm.L[act_i];
                        uint8_t chartowrite = (act_phrase == 0) ? '$' : dict.d[dict.end[act_phrase-1]-w-1];
                        phrase_outlabels[c].push_back(chartowrite);
                        if (tfm.dout[++act_i] == 0) phrase_dout[c].push_back(0);
                        else break;
                    }
                } else {
                    // many incoming edges
                    new_i = tfm.dout_select(din_rank);
                    uint32_t act_phrase = tfm.L[new_i];
                    uint8_t chartowrite = (act_phrase == 0) ? '$' : dict.d[dict.end[act_phrase-1]-w-1];
                    phrase_outlabels[c].push_back(chartowrite);
                    phrase_din[c].push_back(0);
                }
            }
            last_size += phrase_outlabels[c].size();
            bitvector_size += phrase_din[c].size();
            phrase_din[c].push_back(1);
            phrase_dout[c].push_back(1);
        }
        cout << "Phrase outlabels created!" << endl;

        // prepare for serialization
        int_vector<> din (bitvector_size, 0, 1);
        int_vector<> dout (bitvector_size, 0, 1);
        int_vector<> L (last_size, 0, 8);
        uint64_t cnt_din = 0, cnt_dout = 0, cnt_L = 0;
        for (int i = 0; i < dict.dwords+1; i++) {
            cout << "phrase: " << i+1 << endl;
            for (int j = 0; j < phrase_outlabels[i].size(); j++) {
                cout << phrase_outlabels[i][j];
            } cout << endl;
            for (int j = 0; j < phrase_outlabels[i].size(); j++)    L[cnt_L++]          = phrase_outlabels[i][j];
            for (int j = 0; j < phrase_din[i].size(); j++)          din[cnt_din++]      = phrase_din[i][j]; 
            for (int j = 0; j < phrase_dout[i].size(); j++)         dout[cnt_dout++]    = phrase_dout[i][j];
        }

        asprintf(&name, "%s.%s", argv[2], "last");
        store_to_file(L, string (name));
        asprintf(&name, "%s.%s", argv[2], "din");
        store_to_file(din, string(name));
        asprintf(&name, "%s.%s", argv[2], "dout");
        store_to_file(dout, string(name));
        cout << "SOURCES file successfully written!" << endl;
}
