/* **************************************************************************
 * pfbwt.cpp
 * Output the BWT, the SA (option -S) or the sampled SA (option -s)
 * computed using the prefix-free parsing technique
 *  
 * Usage:
 *   pfbwt[NT][64].x -h
 * for usage info
 *
 * See newscan.cpp for a description of what it does  
 * 
 **************************************************************************** */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <semaphore.h>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include "tfm_index.hpp"
extern "C" {
#include "gsa/gsacak.h"
#include "utils.h"
}

using namespace std;
using namespace __gnu_cxx;
using namespace sdsl;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
    char *basename;
    string parseExt =  EXTPARSE;    // extension final parse file  
    string occExt =    EXTOCC;      // extension occurrences file  
    string dictExt =   EXTDICT;     // extension dictionary file  
    string lastExt =   EXTLST;      // extension file containing last chars   
    int w = 10;            // sliding window size and its default 
};

static long get_num_words(uint8_t *d, long n);
static long binsearch(uint_t x, uint_t a[], long n);
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid);
static void compute_dict_bwt_lcp(uint8_t *d, long dsize,long dwords, int w, uint_t **sap, int_t **lcpp);
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, tfm_index<> &tfmp, uint32_t *ilist, FILE *fbwt, long &easy_bwts, long &hard_bwts);

// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
struct SeqId {
    uint32_t id;       // lex. id of the dictionary word to which the suffix belongs
    int remaining;     // remaining copies of the suffix to be considered  
    uint32_t *bwtpos;  // list of bwt positions of this dictionary word
    uint8_t char2write;// char to be written (is the one preceeding the suffix)

    // constructor
    SeqId(uint32_t i, int r, uint32_t *b, int8_t c) : id(i), remaining(r), bwtpos(b) {
        char2write = c;
    }

    // advance to the next bwt position, return false if there are no more positions 
    bool next() {
        remaining--;
        bwtpos += 1;
        return remaining>0;
    }
    bool operator<(const SeqId& a);
};

bool SeqId::operator<(const SeqId& a) {
    return *bwtpos > *(a.bwtpos);
}


inline uint8_t get_prev(int w, uint8_t *d, uint64_t *end, uint32_t seqid) {
    return d[end[seqid]-w-1];
}

void write_bitvector(FILE *f, bool bit, uint8_t &cnt, uint8_t &buffer, bool hard_write=false) {
    buffer |= (bit << (7 - cnt++));
    if (hard_write || (cnt == 8)) {
        if(fputc(buffer,f) == EOF) die("Din/Dout write error 0");
        cnt = 0;
        buffer = 0;
    }
}


/* *******************************************************************
 * Computation of the final BWT
 * 
 * istart[] and islist[] are used together. For each dictionary word i 
 * (considered in lexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i 
 * ******************************************************************* */
void bwt(Args &arg, uint8_t *d, long dsize, uint64_t *end_to_phrase, // dictionary and its size  
         uint32_t *ilist, tfm_index<> &tfmp,long psize, // ilist, last and their size 
         uint32_t *istart, long dwords, uint_t *sa, int_t *lcp) { // starting point in ilist for each word and # words
    
    // set d[0]==0 as this is the EOF char in the final BWT
    assert(d[0]==Dollar);
    d[0] = 0;

    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
    uint_t *eos = sa+1;
    for (int i = 0; i < dwords-1; i++) assert(eos[i] < eos[i+1]);

    // open output file
    FILE *fbwt = open_aux_file(arg.basename,"L","wb");

    // main loop: consider each entry in the SA of dict
    time_t start = time(NULL);
    long full_words = 0, easy_bwts = 0, hard_bwts = 0, next;
    uint32_t seqid;
    for(long i=dwords+arg.w+1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i+1;  // prepare for next iteration  
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        //cout << suffixLen << " " << seqid << endl;
        // ignore suffixes of lenght <= w
        if(suffixLen <= arg.w) continue;
        // ----- simple case: the suffix is a full word 
        if(sa[i] == 0 || d[sa[i]-1] == EndOfWord) {
            full_words++;
            uint32_t start = tfmp.C[seqid+1], end = tfmp.C[seqid+2];
            assert(tfmp.din[start] == 1);
            for (int j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j+1));
                    while (1) {
                        if (tfmp.L[pos] == 0) pos = 0;
                        uint32_t act_phrase = tfmp.L[pos] - 1;
                        uint8_t char_to_write = get_prev(arg.w, d, end_to_phrase, act_phrase);
                        easy_bwts++;
                        //cout << easy_bwts + hard_bwts << " " << seqid << " " << char_to_write << endl;
                        if(fputc(char_to_write,fbwt)==EOF) die("L write error 0");
                        if (tfmp.dout[++pos] == 1) break;
                    }
                }
            }
            continue; // proceed with next i 
        } else {
        // ----- hard case: there can be a group of equal suffixes starting at i
        // save seqid and the corresponding char 
        vector<uint32_t> id2merge(1, seqid); 
        vector<uint8_t> char2write(1,d[sa[i]-1]);
        while(next < dsize && lcp[next] >= suffixLen) {
            assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
            assert(sa[next]>0 && d[sa[next]-1] != EndOfWord); // sa[next] cannot be a full word
            int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
            assert(nextsuffixLen>=suffixLen);
            if(nextsuffixLen==suffixLen) {
                id2merge.push_back(seqid);           // sequence to consider
                char2write.push_back(d[sa[next]-1]);  // corresponding char
                next++;
            }
            else break;
        }
        // output to fbwt the bwt chars corresponding to the current dictionary suffix, and, if requested, some SA values 
        fwrite_chars_same_suffix(id2merge,char2write,tfmp,ilist,fbwt,easy_bwts,hard_bwts);
        }
    }
    assert(full_words == dwords);
    cout << "Full words: " << full_words << endl;
    cout << "Easy bwt chars: " << easy_bwts << endl;
    cout << "Hard bwt chars: " << hard_bwts << endl;
    cout << "Generating the final BWT took " << difftime(time(NULL),start) << " wall clock seconds\n";    
    fclose(fbwt);
}

void din(Args &arg, uint8_t *d, long dsize, uint64_t *end_to_phrase, // dictionary and its size  
         uint32_t *ilist, tfm_index<> &tfmp,long psize, // ilist, last and their size 
         uint32_t *istart, long dwords, uint_t *sa, int_t *lcp) { // starting point in ilist for each word and # words
    
    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
    uint_t *eos = sa+1;
    for (int i = 0; i < dwords-1;i++) assert(eos[i] < eos[i+1]);

    // open output file
    FILE *fdin = open_aux_file(arg.basename,"din","wb");

    // main loop: consider each entry in the SA of dict
    long full_words = 0, easy_bwts = 0, hard_bwts = 0, next;
    uint32_t seqid;
    uint8_t cnt = 0, buffer = 0;
    for(long i=dwords+arg.w+1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i+1;  // prepare for next iteration  
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        //cout << suffixLen << " " << seqid << endl;
        // ignore suffixes of lenght <= w
        if(suffixLen <= arg.w) continue;
        // ----- simple case: the suffix is a full word 
        if(sa[i] == 0 || d[sa[i]-1] == EndOfWord) {
            full_words++;
            uint32_t start = tfmp.C[seqid+1], end = tfmp.C[seqid+2];
            assert(tfmp.din[start] == 1);
            for (int j = start; j < end; j++) {
                write_bitvector(fdin, tfmp.din[j], cnt, buffer);
            }
            continue; // proceed with next i 
        } else {
            // ----- hard case: there can be a group of equal suffixes starting at i
            // save seqid and the corresponding char 
            int bits_to_write = tfmp.C[seqid+2] - tfmp.C[seqid+1];
            while(next < dsize && lcp[next] >= suffixLen) {
                assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
                assert(sa[next]>0 && d[sa[next]-1] != EndOfWord); // sa[next] cannot be a full word
                int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
                assert(nextsuffixLen>=suffixLen);
                if(nextsuffixLen==suffixLen) {
                    bits_to_write += tfmp.C[seqid+2] - tfmp.C[seqid+1];
                    next++;
                }
                else break;
            }
            for (int k = 0; k < bits_to_write; k++) write_bitvector(fdin, 1, cnt, buffer);
        }
    }
    write_bitvector(fdin, 1, cnt, buffer, true);
    fclose(fdin);
}

void dout(Args &arg, uint8_t *d, long dsize, uint64_t *end_to_phrase, // dictionary and its size  
         uint32_t *ilist, tfm_index<> &tfmp,long psize, // ilist, last and their size 
         uint32_t *istart, long dwords, uint_t *sa, int_t *lcp) { // starting point in ilist for each word and # words
    
    // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
    uint_t *eos = sa+1;
    for (int i = 0; i < dwords-1;i++) assert(eos[i] < eos[i+1]);

    // open output file
    FILE *fdout= open_aux_file(arg.basename,"dout","wb");

    // main loop: consider each entry in the SA of dict
    long full_words = 0, easy_bwts = 0, hard_bwts = 0, next;
    uint32_t seqid;
    uint8_t cnt = 0, buffer = 0;
    for(long i=dwords+arg.w+1; i < dsize; i = next) {
        // we are considering d[sa[i]....]
        next = i+1;  // prepare for next iteration  
        // compute length of this suffix and sequence it belongs
        int_t suffixLen = getlen(sa[i], eos, dwords, &seqid);
        //cout << suffixLen << " " << seqid << endl;
        // ignore suffixes of lenght <= w
        if(suffixLen <= arg.w) continue;
        // ----- simple case: the suffix is a full word 
        if(sa[i] == 0 || d[sa[i]-1] == EndOfWord) {
            full_words++;
            uint32_t start = tfmp.C[seqid+1], end = tfmp.C[seqid+2];
            assert(tfmp.din[start] == 1);
            for (int j = start; j < end; j++) {
                if (tfmp.din[j] == 1) {
                    uint32_t pos = tfmp.dout_select(tfmp.din_rank(j+1));
                    if (tfmp.L[pos] == 0) pos = 0;
                    while (1) {
                        write_bitvector(fdout, tfmp.dout[pos], cnt, buffer);
                        if (tfmp.dout[++pos] == 1) break;
                    }
                }
            }
            continue; // proceed with next i 
        } else {
            // ----- hard case: there can be a group of equal suffixes starting at i
            // save seqid and the corresponding char 
            int bits_to_write = tfmp.C[seqid+2] - tfmp.C[seqid+1];
            while(next < dsize && lcp[next] >= suffixLen) {
                assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
                assert(sa[next]>0 && d[sa[next]-1] != EndOfWord); // sa[next] cannot be a full word
                int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
                assert(nextsuffixLen>=suffixLen);
                if(nextsuffixLen==suffixLen) {
                    bits_to_write += tfmp.C[seqid+2] - tfmp.C[seqid+1];
                    next++;
                }
                else break;
            }
            for (int k = 0; k < bits_to_write; k++) write_bitvector(fdout, 1, cnt, buffer);
        }
    }
    write_bitvector(fdout, 1, cnt, buffer, true);
    fclose(fdout);
}

void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-h  \tshow help and exit" << endl;
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
    int c;
    extern char *optarg;
    extern int optind;

    puts("==== Command line:");
    for(int i = 0; i < argc; i++) {
        printf(" %s",argv[i]);
    }
    puts("");

    string sarg;
    while ((c = getopt( argc, argv, "t:w:sehS") ) != -1) {
        switch(c) {
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'h':
            print_help(argv, arg); exit(1);
        case '?':
        cout << "Unknown option. Use -h for help." << endl;
        exit(1);
        }
    }
    // the only input parameter is the file name
    arg.basename = NULL; 
    if (argc == optind+1) 
        arg.basename = argv[optind];
    else {
        cout << "Invalid number of arguments" << endl;
        print_help(argv,arg);
    }
    if(arg.w <2) {
        cout << "Windows size must be at least 2\n";
        exit(1);
    }
}


int main(int argc, char** argv) {
    time_t start = time(NULL);  

    // translate command line parameters
    Args arg;
    parseArgs(argc, argv, arg);
    // read dictionary file 
    FILE *g = open_aux_file(arg.basename,EXTDICT,"rb");
    fseek(g,0,SEEK_END);
    long dsize = ftell(g);
    if(dsize<0) die("ftell dictionary");
    if(dsize<=1+arg.w) die("invalid dictionary file");
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
    long e = fread(d,1,dsize,g);
    if(e!=dsize) die("fread");
    fclose(g);

    uint32_t dwords = get_num_words(d, dsize);
    
    // map end of phrases to phrase ids
    uint64_t *end_to_phrase = new uint64_t[dwords], cnt = 0;
    for (int i = 0; i < dsize; i++) {
        if (d[i] == EndOfWord) end_to_phrase[cnt++] = i;
    }

    // read ilist file 
    g = open_aux_file(arg.basename,EXTILIST,"rb");
    fseek(g,0,SEEK_END);
    e = ftell(g);
    if(e<0) die("ftell ilist file");
    if(e%4!=0) die("invalid ilist file");
    long psize = e/4;
    cout  << "Parsing size: " << psize << endl;
    if(psize>0xFFFFFFFEL) die("More than 2^32 -2 words in the parsing");
    uint32_t *ilist = new uint32_t[psize];  
    rewind(g);
    e = fread(ilist,4,psize,g);
    if(e!=psize) die("fread 3");
    fclose(g);

    // read the tunneled WG of the parse
    char *name;
    asprintf(&name, "%s.%s", arg.basename, "tunnel");
    tfm_index<> tfmp;
    load_from_file(tfmp, name);

    // convert occ entries into starting positions inside ilist
    // ilist also contains the position of EOF but we don't care about it since it is not in dict 
    uint32_t *occ = new uint32_t[dwords+1];
    uint32_t last=0; // starting position in ilist of the smallest dictionary word  
    for(long i=0;i<dwords;i++) {
        uint32_t tmp = tfmp.C[i+2] - tfmp.C[i+1];
        occ[i] = last;
        last += tmp;
        //cout << i << " = " << occ[i] << " " << tfmp.C[i+1] - tfmp.C[i] << endl; 
    }
    occ[dwords]=last;

    // compute SA and BWT of D and do some checking on them 
    uint_t *sa; int_t *lcp;
    compute_dict_bwt_lcp(d, dsize, dwords, arg.w, &sa, &lcp);

    bwt(arg,d,dsize,end_to_phrase,ilist,tfmp,psize,occ,dwords,sa,lcp);
    din(arg,d,dsize,end_to_phrase,ilist,tfmp,psize,occ,dwords,sa,lcp);
    dout(arg,d,dsize,end_to_phrase,ilist,tfmp,psize,occ,dwords,sa,lcp);
    delete[] ilist;
    delete[] occ;
    delete[] d;  
    delete[] lcp;
    delete[] sa;
    cout << "==== Elapsed time: " << difftime(time(NULL),start) << " wall clock seconds\n";      
    return 0;
}

// --------------------- aux functions ----------------------------------

// compute the number of words in a dictionary
static long get_num_words(uint8_t *d, long n) {
    long i, num = 0;
    for(i = 0; i < n; i++) if(d[i]==EndOfWord) num++;
    assert(d[n-1]==EndOfDict);
    return num;
}

// binary search for x in an array a[0..n-1] that doesn't contain x
// return the lowest position that is larger than x
static long binsearch(uint_t x, uint_t a[], long n) {
    long lo=0; long hi = n-1;
    while(hi>lo) {
        assert( ((lo==0) || x>a[lo-1]) && x< a[hi]);
        int mid = (lo+hi)/2;
        assert(x!=a[mid]);  // x is not in a[]
        if(x<a[mid]) hi = mid;
        else lo = mid+1;
    }
    assert(((hi==0) || x>a[hi-1]) && x< a[hi]);
    return hi; 
}


// return the length of the suffix starting in position p.
// also write to seqid the id of the sequence containing that suffix 
// n is the # of distinct words in the dictionary, hence the length of eos[]
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid) {
    assert(p<eos[n-1]);
    *seqid = binsearch(p,eos,n);
    assert(eos[*seqid]> p); // distance between position p and the next $
    return eos[*seqid] - p;
}

// compute the SA and LCP array for the set of (unique) dictionary words
// using gSACA-K. Also do some checking based on the number and order of the special symbols
// d[0..dsize-1] is the dictionary consisting of the concatenation of dictionary words
// in lex order with EndOfWord (0x1) at the end of each word and 
// d[size-1] = EndOfDict (0x0) at the very end. It is d[0]=Dollar (0x2)
// since the first words starts with $. There is another word somewhere
// ending with Dollar^wEndOfWord (it is the last word in the parsing,
// but its lex rank is unknown).  
static void compute_dict_bwt_lcp(uint8_t *d, long dsize,long dwords, int w, 
                          uint_t **sap, int_t **lcpp) // output parameters
{
    uint_t *sa = new uint_t[dsize];
    int_t *lcp = new int_t[dsize];
    (void) dwords; (void) w;

    cout  << "Each SA entry: " << sizeof(*sa) << " bytes\n";
    cout  << "Each LCP entry: " << sizeof(*lcp) << " bytes\n";

    cout << "Computing SA and LCP of dictionary" << endl; 
    time_t  start = time(NULL);
    gsacak(d,sa,lcp,NULL,dsize);
    cout << "Computing SA/LCP took " << difftime(time(NULL),start) << " wall clock seconds\n";  
    // ------ do some checking on the sa
    assert(d[dsize-1]==EndOfDict);
    assert(sa[0]==(unsigned long)dsize-1);// sa[0] is the EndOfDict symbol 
    for(long i=0;i<dwords;i++) assert(d[sa[i+1]]==EndOfWord); // there are dwords EndOfWord symbols 
    // EndOfWord symbols are in position order, so the last is d[dsize-2]    
    assert(sa[dwords]==(unsigned long)dsize-2);  
    // there are wsize+1 $ symbols: 
    // one at the beginning of the first word, wsize at the end of the last word
    for(long i=0;i<=w;i++) assert(d[sa[i+dwords+1]]==Dollar);         
    // in sa[dwords+w+1] we have the first word in the parsing since that $ is the lex. larger  
    assert(d[0]==Dollar);
    assert(sa[dwords+w+1]==0);
    assert(d[sa[dwords+w+2]]>Dollar);  // end of Dollar chars in the first column
    assert(lcp[dwords+w+2]==0); 
    // copy sa and lcp address
    *sap = sa;  *lcpp = lcp;  
}


// write to the bwt all the characters preceding a given suffix
// doing a merge operation if necessary
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, tfm_index<> &tfmp,
                                    uint32_t *ilist, FILE *fbwt, long &easy_bwts, long &hard_bwts) {
    size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix
    bool samechar = true;
    for(size_t i = 1; (i < numwords) && samechar; i++) {
        samechar = (char2write[i-1]==char2write[i]);
    }

    if(samechar) {
        for(size_t i = 0; i < numwords; i++) {
            uint32_t s = id2merge[i]+1;
            for(long j = tfmp.C[s]; j < tfmp.C[s+1]; j++) {
                if(fputc(char2write[0],fbwt) == EOF) die("L write error 1");
            }
            easy_bwts +=  tfmp.C[s+1] - tfmp.C[s]; 
            //cout << easy_bwts + hard_bwts  << ":: numwords = " << numwords  <<" " << char2write[0] << " " << tfmp.C[s+1] - tfmp.C[s] << endl;
        }
    } else {  // many words, many chars...     
        //cout << "STARTING SORTING " << endl;
        vector<SeqId> heap; // create heap
        for(size_t i=0; i<numwords; i++) {
            uint32_t s = id2merge[i]+1;
            //cout << "phrase: " << s << " pos: " << ilist[tfmp.C[s]-1] << endl;
            heap.push_back(SeqId(s,tfmp.C[s+1]-tfmp.C[s], ilist+(tfmp.C[s]-1), char2write[i]));
        }
        std::make_heap(heap.begin(),heap.end());
        while(heap.size()>0) {
            // output char for the top of the heap
            SeqId s = heap.front();
            if(fputc(s.char2write,fbwt)==EOF) die("L write error 2");
            hard_bwts += 1;
            // remove top 
            //cout << "HEAPSORT:" <<  easy_bwts + hard_bwts  << ":::" << (*s.bwtpos) << "--" << s.char2write << " " << endl;
            pop_heap(heap.begin(),heap.end());
            heap.pop_back();
            // if remaining positions, reinsert to heap
            if(s.next()) {
                heap.push_back(s);
                push_heap(heap.begin(),heap.end());
            }
        }
    }
}

