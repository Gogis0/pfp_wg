/* ******************************************************************************
 * newscan.cpp
 * 
 * parsing algorithm for bwt construction of repetitive sequences based 
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/)
 * 
 * Usage:
 *   newscan.x wsize modulus file
 * 
 * Unless the parameter -c (compression rather than BWT construction,
 * see "Compression mode" below) is given, the input file cannot contain 
 * the characters 0x0, 0x1, 0x2 which are used internally for BWT construction 
 * 
 * Since the i-th thread accesses the i-th segment of the input file 
 * random access (fseek) must be possible. For gzipped inputs use 
 * cnewscan.x which doesn't use threads but automatically extracts the 
 * content from a gzipped input using the lz library. 
 * 
 * The parameters wsize and modulus are used to define the prefix free parsing 
 * using KR-fingerprints (see paper)
 * 
 * 
 * *** BWT construction ***
 *  
 * The algorithm computes the prefix free parsing of 
 *     T = (0x2)file_content(0x2)^wsize
 * in a dictionary of words D and a parsing P of T in terms of the  
 * dictionary words. Note that the words in the parsing overlap by wsize.
 * Let d denote the number of words in D and p the number of phrases in 
 * the parsing P
 * 
 * newscan.x outputs the following files:
 * 
 * file.dict
 * containing the dictionary words in lexicographic order with a 0x1 at the end of
 * each word and a 0x0 at the end of the file. Size: |D| + d + 1 where
 * |D| is the sum of the word lengths
 * 
 * file.occ
 * the number of occurrences of each word in lexicographic order.
 * We assume the number of occurrences of each word is at most 2^32-1
 * so the size is 4d bytes
 * 
 * file.parse
 * containing the parse P with each word identified with its 1-based lexicographic 
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes
 * 
 * file.last 
 * containing the character in position w+1 from the end for each dictionary word
 * Size: p
 * 
 * file.sai (if option -s is given on the command line) 
 * containing the ending position +1 of each parsed word in the original
 * text written using IBYTES bytes for each entry (IBYTES defined in utils.h) 
 * Size: p*IBYTES
 * 
 * The output of newscan.x must be processed by bwtparse, which invoked as
 * 
 *    bwtparse file
 * 
 * computes the BWT of file.parse and produces file.ilist of size 4p+4 bytes
 * contaning, for each dictionary word in lexicographic order, the list 
 * of BWT positions where that word appears (ie i\in ilist(w) <=> BWT[i]=w).
 * There is also an entry for the EOF word which is not in the dictionary 
 * but is assumed to be the smallest word.  
 * 
 * In addition, bwtparse permutes file.last according to
 * the BWT permutation and generates file.bwlast such that file.bwlast[i] 
 * is the char from P[SA[i]-2] (if SA[i]==0 , BWT[i]=0 and file.bwlast[i]=0, 
 * if SA[i]==1, BWT[i]=P[0] and file.bwlast[i] is taken from P[n-1], the last 
 * word in the parsing).  
 * 
 * If the option -s is given to bwtparse, it permutes file.sai according
 * to the BWT permutation and generate file.bwsai using again IBYTES
 * per entry.  file.bwsai[i] is the ending position+1 of BWT[i] in the 
 * original text 
 * 
 * The output of bwtparse (the files .ilist .bwlast) together with the
 * dictionary itself (file .dict) and the number of occurrences
 * of each word (file .occ) are used to compute the final BWT by the 
 * pfbwt algorithm.
 *
 * As an additional check to the correctness of the parsing, it is 
 * possible to reconstruct the original file from the files .dict
 * and .parse using the unparse tool.  
 * 
 * 
 *  *** Compression mode ***
 * 
 * If the -c option is used, the parsing is computed for compression
 * purposes rather than for building the BWT. In this case the redundant 
 * information (phrases overlaps and 0x2's) is not written to the output files.
 * 
 * In addition, the input can contain also the characters 0x0, 0x1, 0x2
 * (ie can be any input file). The program computes a quasi prefix-free 
 * parsing (with no overlaps): 
 * 
 *   T = w_0 w_1 w_2 ... w_{p-1}
 *
 * where each word w_i, except the last one, ends with a lenght-w suffix s_i
 * such that KR(s_i) mod p = 0 and s_i is the only lenght-w substring of
 * w_i with that property, with the possible exception of the lenght-w
 * prefix of w_0.
 * 
 * In Compression mode newscan.x outputs the following files:
 * 
 * file.dicz
 * containing the concatenation of the (distinct) dictionary words in 
 * lexicographic order. 
 * Size: |D| where |D| is the sum of the word lengths
 * 
 * file.dicz.len
 * containing the lenght in bytes of the dictionary words again in 
 * lexicographic order. Each lenght is represented by a 32 bit int.
 * Size: 4d where d is the number of distinct dictionary words. 
 * 
 * file.parse
 * containing the parse P with each word identified with its 1-based lexicographic 
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes.
 * 
 * From the above three files it is possible to recover the original input
 * using the unparsz tool.
 * 
 * 
 * *** Usage of alphabets larger than 256 ***
 * 
 * The -b B command line option -b makes it possible to use symbols consisting of B bytes
 * The size of the input file must be a multiple of B. The algorithm will produce 
 * a dictionary consisting of words made of B-byte symbol; the dictionary will still be
 * lexicographically sorted. The program will assume that the bytes inside each symbol
 * are in Little Endian format, unless the -g command line option is used.  
 * Note that the length of the input and of the dictionary words will be 
 * expressed in symbols.
 * The -b option is supported for both compression mode and BWT construction
 * but at the moment the other stages of BWT construction do not support 
 * alphabets larger than 256
 */
#include "newscan_functions.cpp"

uint64_t my_process_file(Args &arg, map<uint64_t, word_stats> &wordFreq) {

    ifstream f(arg.inputFileName);
    if (!f.rdbuf()->is_open()) {
        perror("can not open file for reading");
        exit(EXIT_FAILURE);
    }

    FILE *g = open_aux_file(arg.inputFileName.c_str(), "parse_old", "wb");
    // FILE *last_file = open_aux_file(arg.inputFileName.c_str(), "last", "wb");
    uint8_t buffer[arg.bytexsymb];
    uint8_t dollar[arg.bytexsymb];
    dollar[arg.bytexsymb - 1] = Dollar;  // this is the generalized Dollar symbol

    // main loop on the symbols of the input file
    uint64_t pos = 0; // ending position +1 of previous word in the original text, used for computing sa_info
    assert(IBYTES <= sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
    KR_window krw(arg.w, arg.bytexsymb);
    ztring word;
    ztring_append(word, dollar, arg.bytexsymb);

    while (true) {
        f.read((char *) buffer, arg.bytexsymb);
        // we must be able to read exactly arg.bytexsymb bytes otherwise the symbol is incomplete
        if (f.gcount() == 0) break;
        if (f.gcount() != arg.bytexsymb) {
            cerr << "Incomplete symbol at position " << f.tellg() - f.gcount() << " Exiting....\n";
            exit(1);
        }
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
            save_update_word(arg, word, wordFreq, g, NULL, NULL, pos);
        }
    }
    // check that we really reached the end of the file
    if (!f.eof()) die("Error reading from input file (process_file)");
    // virtually add w null chars at the end of the file and add the last word in the dict
    for (int i = 0; i < arg.w; i++) ztring_append(word, dollar, arg.bytexsymb);
    save_update_word(arg, word, wordFreq, g, NULL, NULL, pos);

    // close input and output files
    // if (last_file) if (fclose(last_file) != 0) die("Error closing last file");
    if (fclose(g) != 0) die("Error closing parse file");
    assert(pos == krw.tot_symb + arg.w);
    f.close();
    return krw.tot_symb;
}

int main() {
    Args arg;
    arg.inputFileName = "../data/yeast.fasta";
    arg.w = 4;
    arg.p = 11;
    map<uint64_t, word_stats> wordFreq;

    my_process_file(arg, wordFreq);

    uint64_t totDWord = wordFreq.size();
    vector<const ztring *> dictArray;
    dictArray.reserve(totDWord);
    uint64_t sumLen = 0;
    uint64_t totWord = 0;
    for (auto &x: wordFreq) {
        sumLen += (x.second.str.size()) / arg.bytexsymb;
        totWord += x.second.occ;
        dictArray.push_back(&x.second.str);
    }
    assert(dictArray.size() == totDWord);
    sort(dictArray.begin(), dictArray.end(), pztringCompare);
    writeDictOcc(arg, wordFreq, dictArray);
    dictArray.clear();
    remapParse(arg, wordFreq);
    return 0;
}