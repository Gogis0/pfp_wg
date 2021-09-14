#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>

extern "C" {
#include "../bigrepair/ctph/utils.h"
}

using namespace std;
using namespace __gnu_cxx;



// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;
// type used to represent chars when PARSE_WORDS is defined
#ifdef PARSE_WORDS
typedef uint32_t char_int_t;
#endif
// type used to represent a sequence of input symbol trasformed to bytes
typedef vector<uint8_t> ztring;


// values of the wordFreq map: word, its number of occurrences, and its lex rank
struct word_stats {
    ztring str;
    occ_int_t occ;
    word_int_t rank = 0;
};


// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
    string inputFileName;
    int w = 10;            // sliding window size and its default
    int p = 100;           // modulus for establishing stopping w-tuples
    int bytexsymb = 1;     // number of bytes per symbol
    bool bigEndian = false;// when bytexsymb>1 whether they are bigEndian or smallEndian
    bool SAinfo = false;   // compute SA information
    bool compress = false; // parsing called in compress mode
    int th = 0;              // number of helper threads
    int verbose = 0;         // verbosity level
};


// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
    int wsize;            // number of symbols in window
    int bytexsymb;        // number of bytes per symbol
    int wbsize;           // size of window in bytes
    uint8_t *window;
    int asize;            // alphabet size
    const uint64_t prime = 1999999973; // slightly less that 2^31
    uint64_t hash;        // hash of the symbols currently in window
    uint64_t tot_symb;    // char added so far, equals symbols*bxs
    uint64_t asize_pot;   // asize^(wsize-1) mod prime

    KR_window(int w, int bxs) : wsize(w), bytexsymb(bxs) {
        wbsize = wsize * bytexsymb;    // size of window in bytes
        asize = 256;                 // alphabet size for bytes
        asize_pot = modpow(asize, wbsize - 1); // power used to update hash when oldest char exit
        // alloc and clear window
        window = new uint8_t[wbsize];
        reset();
    }

    // power modulo prime
    uint64_t modpow(uint64_t base, uint64_t exp) {
        assert(exp > 0);
        if (exp == 1) return base;
        if (exp % 2 == 0)
            return modpow((base * base) % prime, exp / 2);
        else
            return (base * modpow((base * base) % prime, exp / 2)) % prime;
    }

    // init window, hash, and tot_symb
    void reset() {
        for (int i = 0; i < wbsize; i++) window[i] = 0;
        // init hash value and related values
        hash = tot_symb = 0;
    }

    // add a symbol consisting of bytexsymbol uint8's to the window
    // update and return the hash for the resulting window
    uint64_t addsymbol(uint8_t *s) {
        // compute destination of symbol's bytes inside window[]
        int k = (tot_symb++ % wsize) * bytexsymb;
        assert(k + bytexsymb - 1 < wbsize); // make sure we are inside window[]
        for (int i = 0; i < bytexsymb; i++) {
            // complex expression to avoid negative numbers
            hash += (prime - (window[k] * asize_pot) % prime); // remove window[k] contribution
            hash = (asize * hash + s[i]) % prime;      //  add char i
            window[k++] = s[i];
        }
        return hash;
    }

    ~KR_window() {
        delete[] window;
    }

};
// -----------------------------------------------------------

// ---- append n bytes of b[] to the ztring w; used to add a new symbol to the current word
static void ztring_append(ztring &w, uint8_t *b, int n);

// ---- reverse a byte buffer
static void buffer_reverse(uint8_t *b, int n);

// ---- transform integers to littelEndian anb write them to file
void fwrite_littleEndian(const uint8_t *b, int s, int n, FILE *f);

// ---- process a compete word
static void
save_update_word(Args &arg, ztring &w, map<uint64_t, word_stats> &freq, FILE *tmp_parse_file, FILE *last, FILE *sa,
                 uint64_t &pos);

// ----------- overloading << for vectors
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
        out << '[';
        for (auto i: v) out << i << " ";
        out << "\b]";
    }
    return out;
}

// compute 64-bit KR hash of a ztring
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(ztring s) {
    uint64_t hash = 1;    // forcing a nonzero MSB to avoid collisions between 000xyz and xyz
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for (size_t k = 0; k < s.size(); k++) {
        int c = s[k];
        assert(c >= 0 && c < 256);
        hash = (256 * hash + c) % prime;    //  add char k
    }
    return hash;
}


// save current word w in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void
save_update_word(
    Args &arg, ztring &w, map<uint64_t, word_stats> &freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos
){
    size_t minsize = arg.w;  // no word should be smaller than windows size
    assert(w.size() % arg.bytexsymb == 0);
    assert(pos == 0 || w.size() > (minsize * arg.bytexsymb));
    if (w.size() <= minsize * arg.bytexsymb) return;
    // save overlap consisting of the last minsize symbols
    ztring overlap(w.begin() + (w.size() - minsize * arg.bytexsymb), w.end());
    // if we are compressing, discard the overlap (the last minsize symbols) from current word before storage

    // get the hash value and write it to the temporary parse file
    uint64_t hash = kr_hash(w);
    if (fwrite(&hash, sizeof(hash), 1, tmp_parse_file) != 1) die("parse write error");

    // update frequency table for current hash
    if (freq.find(hash) == freq.end()) {
        freq[hash].occ = 1; // new hash
        freq[hash].str = w;
    } else {
        freq[hash].occ += 1; // known hash
        if (freq[hash].occ <= 0) {
            cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
            cerr << MAX_WORD_OCC << ") exceeded\n";
            exit(1);
        }
        if (freq[hash].str != w) {
            cerr << "Emergency exit! Hash collision for strings:\n";
            cerr << freq[hash].str << "\n  vs\n" << w << endl;
            exit(1);
        }
    }

    if (arg.compress)
        pos += w.size() / arg.bytexsymb; // if compressing, just update position
    // keep only the overlapping part of the window
    w.assign(overlap.begin(), overlap.end());
    assert((int) w.size() == arg.w * arg.bytexsymb);
}

// prefix free parse of a file. the main input parameters are in arg
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args &arg, map<uint64_t, word_stats> &wordFreq) {
    //open a, possibly compressed, input file
    string fnam = arg.inputFileName;
    ifstream f(fnam);
    if (!f.rdbuf()->is_open()) {// is_open does not work on igzstreams
        perror(__func__);
        throw new std::runtime_error("Cannot open input file " + fnam);
    }

    // open the 1st pass parsing file
    FILE *g = open_aux_file(arg.inputFileName.c_str(), EXTPARS0, "wb");
    FILE *sa_file = NULL, *last_file = NULL;
    if (!arg.compress) {
        // open output file containing the symbol at position -(w+1) of each word
        last_file = open_aux_file(arg.inputFileName.c_str(), EXTLST, "wb");
        // if requested open file containing the ending position+1 of each word
        if (arg.SAinfo)
            sa_file = open_aux_file(arg.inputFileName.c_str(), EXTSAI, "wb");
    }

    // init buffers containing a single symbol
    assert(arg.bytexsymb > 0);
    uint8_t buffer[arg.bytexsymb], dollar[arg.bytexsymb] = {0};
    dollar[arg.bytexsymb - 1] = Dollar;  // this is the generalized Dollar symbol

    // main loop on the symbols of the input file
    uint64_t pos = 0; // ending position +1 of previous word in the original text, used for computing sa_info
    assert(IBYTES <= sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
    // init empty KR window
    KR_window krw(arg.w, arg.bytexsymb);
    // init first word in the parsing with a generalized dollar symbol unless we are just compressing
    ztring word;
    if (!arg.compress)
        ztring_append(word, dollar, arg.bytexsymb);

    while (true) {
        f.read((char *) buffer, arg.bytexsymb);
        // we must be able to read exactly arg.bytexsymb bytes otherwise the symbol is incomplete
        if (f.gcount() != arg.bytexsymb) {
            if (f.gcount() == 0) break;
            else {
                cerr << "Incomplete symbol at position " << f.tellg() - f.gcount() << " Exiting....\n";
                exit(1);
            }
        }
        // if not bigEndian swap bytes as we compare byte sequences
        if (arg.bytexsymb > 1 && !arg.bigEndian)
            buffer_reverse(buffer, arg.bytexsymb);
        // if we are not simply compressing then we cannot accept 0,1,or 2
        if (!arg.compress && memcmp(buffer, dollar, arg.bytexsymb) <= 0) {
            cerr << "Invalid symbol (starting with " << buffer[0] << ") at file position ";
            cerr << ((long) f.tellg() - arg.bytexsymb) << ". Exiting...\n";
            exit(1);
        }
        // add new symbol to current word and check if we reached a splitting point
        ztring_append(word, buffer, arg.bytexsymb);
        uint64_t hash = krw.addsymbol(buffer);
        if (hash % arg.p == 0) {
            // end of word, save it and write its full hash to the output file
            save_update_word(arg, word, wordFreq, g, last_file, sa_file, pos);
        }
    }
    // check that we really reached the end of the file
    if (!f.eof()) die("Error reading from input file (process_file)");
    // virtually add w null chars at the end of the file and add the last word in the dict
    for (int i = 0; i < arg.w; i++)
        ztring_append(word, dollar, arg.bytexsymb);
    save_update_word(arg, word, wordFreq, g, last_file, sa_file, pos);

    // close input and output files
    if (sa_file) if (fclose(sa_file) != 0) die("Error closing SA file");
    if (last_file) if (fclose(last_file) != 0) die("Error closing last file");
    if (fclose(g) != 0) die("Error closing parse file");
    if (arg.compress) {
        if (pos != krw.tot_symb) cerr << "pos: " << f.tellg() << " tot_symb " << krw.tot_symb << endl;
        assert(pos == krw.tot_symb);
    } else
        assert(pos == krw.tot_symb + arg.w);
    f.close();
    return krw.tot_symb;
}


// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t, word_stats> &wfreq, vector<const ztring *> &sortedDict) {
    assert(sortedDict.size() == wfreq.size());
    FILE *fdict, *fwlen = NULL, *focc = NULL;
    // open dictionary and occ files
    if (arg.compress) {
        fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICZ, "wb");
        fwlen = open_aux_file(arg.inputFileName.c_str(), EXTDZLEN, "wb");
    } else {
        fdict = open_aux_file(arg.inputFileName.c_str(), EXTDICT, "wb");
        focc = open_aux_file(arg.inputFileName.c_str(), EXTOCC, "wb");
    }

    word_int_t wrank = 1; // current word rank (1 based)
    for (auto x: sortedDict) {                 // *x is the ztring representing the dictionary word
        const uint8_t *word = (*x).data();      // byte sequence of current dictionary word
        assert((*x).size() > 0 && (*x).size() % arg.bytexsymb == 0);
        size_t len = (*x).size() / arg.bytexsymb; // length of word in symbols
        assert(len > (size_t) arg.w || arg.compress);
        uint64_t hash = kr_hash(*x);
        auto &wf = wfreq.at(hash);
        assert(wf.occ > 0);
        //  write word to dictionary
        if (arg.bytexsymb > 1 && !arg.bigEndian) {
            // necessary to transform back the integers to littleEndian
            fwrite_littleEndian(word, arg.bytexsymb, len, fdict);
        } else { // straight copy
            size_t s = fwrite(word, arg.bytexsymb, len, fdict);
            if (s != len) die("Error writing to DICT file");
        }
        // if in compress mode write len as an int32 in a separate file
        if (arg.compress) {
            size_t s = fwrite(&len, 4, 1, fwlen);
            if (s != 1) die("Error writing to WLEN file");
        }
            // otherwise write a EndOfWord symbol on the same file and update occ file
        else {
            // create EndOfWord symbol
            uint8_t tmp[arg.bytexsymb] = {0};
            tmp[arg.bigEndian ? arg.bytexsymb - 1 : 0] = EndOfWord;
            int s = fwrite(tmp, 1, arg.bytexsymb, fdict); // append EndOfWord to the dictionary
            if (s != arg.bytexsymb) die("Error writing EndOfWord to DICT file");
            // write number of occ of current word to focc file
            s = fwrite(&wf.occ, sizeof(wf.occ), 1, focc);
            if (s != 1) die("Error writing to OCC file");
        }
        assert(wf.rank == 0);
        wf.rank = wrank++;
    }
    if (arg.compress) {
        if (fclose(fwlen) != 0) die("Error closing WLEN file");
    } else {
        if (fputc(EndOfDict, fdict) == EOF) die("Error writing EndOfDict to DICT file");
        if (fclose(focc) != 0) die("Error closing OCC file");
    }
    if (fclose(fdict) != 0) die("Error closing DICT file");
}


// remap the parse file replacing hash with rank in the sorted dictionary
void remapParse(Args &arg, map<uint64_t, word_stats> &wfreq) {
    // open parse files. the old parse has been stored in a single file or in multiple files
    mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
    FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");

    // recompute occ as an extra check
    vector<occ_int_t> occ(wfreq.size() + 1, 0); // ranks are zero based
    uint64_t hash;
    while (true) {
        size_t s = mfread(&hash, sizeof(hash), 1, moldp);
        if (s == 0) break;
        if (s != 1) die("Unexpected parse EOF");
        word_int_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        s = fwrite(&rank, sizeof(rank), 1, newp);
        if (s != 1) die("Error writing to new parse file");
    }
    if (fclose(newp) != 0) die("Error closing new parse file");
    if (mfclose(moldp) != 0) die("Error closing old parse segment");
    // check old and recomputed occ coincide
    for (auto &x : wfreq)
        assert(x.second.occ == occ[x.second.rank]);
}


void print_help(char **argv, Args &args) {
    cout << "Usage: " << argv[0] << " <input filename> [options]" << endl;
    cout << "  Options: " << endl
         << "\t-w W\tsliding window size, def. " << args.w << endl
         << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
         << "\t-b B\tnumber of bytes x symbol, def. " << args.bytexsymb << endl
         << "\t-g  \tsymbols are in bigEndian format, def No" << endl
         << "\t-c  \tcompression mode: discard overlaps and $'s" << endl
         << "\t-h  \tshow help and exit" << endl
         << "\t-s  \tcompute suffix array info" << endl;
    exit(1);
}

void parseArgs(int argc, char **argv, Args &arg) {
    int c;
    string sarg;
    while ((c = getopt(argc, argv, "p:w:b:gsht:vc")) != -1) {
        switch (c) {
            case 's':
                arg.SAinfo = true;
                break;
            case 'c':
                arg.compress = true;
                break;
            case 'w':
                sarg.assign(optarg);
                arg.w = stoi(sarg);
                break;
            case 'p':
                sarg.assign(optarg);
                arg.p = stoi(sarg);
                break;
            case 'b':
                sarg.assign(optarg);
                arg.bytexsymb = stoi(sarg);
                break;
            case 'g':
                arg.bigEndian = true;
                break;
            case 't':
                sarg.assign(optarg);
                arg.th = stoi(sarg);
                break;
            case 'v':
                arg.verbose++;
                break;
            case 'h':
                print_help(argv, arg);
                exit(1);
            case '?':
                cout << "Unknown option. Use -h for help." << endl;
                exit(1);
        }
    }
    // the only mandatory input parameter is the file name
    if (argc == optind + 1) {
        arg.inputFileName.assign(argv[optind]);
    } else {
        cout << "Invalid number of arguments" << endl;
        print_help(argv, arg);
    }
    // check algorithm parameters
    if (arg.w < 4) {
        cout << "Windows size must be at least 4\n";    // why?
        exit(1);
    }
    if (arg.p < 10) {
        cout << "Modulus must be at least 10\n";        // why?
        exit(1);
    }
    if (arg.bytexsymb < 1) {
        cout << "# bytes x symbol must be positive\n";
        exit(1);
    }
    if (arg.bytexsymb == 1 && arg.bigEndian) {
        cout << "Option big Endian makes sense only if byte x symbol > 1\n";
        exit(1);
    }
    if (arg.th != 0) {
        cout << "The NT version cannot use threads\n";
        exit(1);
    }
}


// function used to compare two ztring pointers, used to lex sort dictionary words
static bool pztringCompare(const ztring *a, const ztring *b) {
    return *a <= *b;
}

// ---- reverse a byte buffer
static void buffer_reverse(uint8_t *b, int n) {
    int i, j;
    for (i = 0, j = n - 1; i < j; i++, j--) {
        uint8_t t = b[i];
        b[i] = b[j];
        b[j] = t;
    }
}


// write a block of n symbols each one taking size bytes to file f
// first reversing the bytes of each symbol
void fwrite_littleEndian(const uint8_t *b, int size, int n, FILE *f) {
    assert(size > 1 && n > 0);
    uint8_t *tmp = new uint8_t[n * size];

    // each symbol takes a contiguous block of bytes
    // with remainder 0,1,...size-1
    // we remap them to remainder size-1,...1,0
    for (int i = 0; i < n * size; i++)
        tmp[(i / size) * size + (size - 1 - (i % size))] = b[i];

    int s = fwrite(tmp, size, n, f);
    if (s != n) die("Error writing to DICT file (fwrite_littleEndian");
    delete[] tmp;
}


// append n bytes of b[] to the ztring w
// used to add a new symbol to the current word
static void ztring_append(ztring &w, uint8_t *b, int n) {
    for (int i = 0; i < n; i++)
        w.push_back(b[i]);
}


