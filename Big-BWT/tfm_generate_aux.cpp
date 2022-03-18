/*
 * tfm_index_invert.cpp for Edge minimization in de Bruijn graphs
 * Copyright (c) 2019 Uwe Baier, Pascal Weber All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
	cerr << "USAGE: " << argv[0] << "w TFMFILE OCCFILE" << endl;
	cerr << "TFMFILE:" << endl;
	cerr << "  File where to store the serialized trie" << endl;
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

int main(int argc, char **argv) {
	//check parameters
	if (argc < 1) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}
        int w = atoi(argv[5]);
        struct Dict dict = read_dictionary(argv[1]);

	//load tunneled fm index
	tfm_index<> tfm;
	load_from_file( tfm, argv[2] );
        vector<vector<uint8_t>> phrase_sources(tfm.dout_rank(tfm.L.size())); // the list of ourgoing edges for each vertex
        vector<int> phrase_occs(dict.dwords, 0);

        int node = 0;
        uint32_t c = -1;
	for (size_type i = 0; i < tfm.L.size()-1; i++) {
            c = tfm.L[i];
            phrase_occs[c]++;
            //cout << dict.d[1] << endl;
            if (c == 0) phrase_sources[node].push_back('$');
            else phrase_sources[node].push_back(dict.d[dict.end[c-1]-w-1]);
            if (tfm.dout[i+1] == 1) node++;
	}
        phrase_sources[node].push_back(dict.d[dict.end[c-1]-w-1]);
        cout << "search completed" << endl;

        // write the numbers of occurrences for each phrase
        FILE *focc = fopen(argv[3], "wb");
        for (int i = 0; i < dict.dwords; i++) {
            size_t s = fwrite(&phrase_occs[i], sizeof(phrase_occs[i]), 1, focc);
            if (s != 1) {
                cout << "Error writing to OCC file" << endl;
                exit(1);
            }
        }
        cout << "OCC file successfully written!" << endl;

        // write the numbers of occurrences for each phrase
        FILE *fsrc = fopen(argv[4], "wb");
        for (int i = 0; i < phrase_sources.size(); i++) {
            cout << i << " : ";
            for (int j = 0; j < phrase_sources[i].size(); j++) {
                cout << (char)phrase_sources[i][j] << ' ';
                size_t s = fwrite(&phrase_sources[i], sizeof(phrase_sources[i]), 1, fsrc);
                if (s != 1) {
                    cout << "Error writing to SOURCES file" << endl;
                    exit(1);
                }
            }
            cout << endl;
        }
        cout << "SOURCES file successfully written!" << endl;

        delete dict.d;
}
