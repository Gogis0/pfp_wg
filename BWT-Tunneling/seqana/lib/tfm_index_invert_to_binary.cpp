#include <iostream>
#include <deque>
#include "../include/tfm_index.hpp"

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage(char **argv) {
    cerr << "USAGE: " << argv[0] << " <tunneled parse> <parse>" << endl;
    cerr << "TFMFILE:" << endl;
    cerr << "  File where to store the serialized trie" << endl;
    exit(EXIT_FAILURE);
};

int main(int argc, char **argv) {
    if (argc != 3) printUsage(argv);
    //load tunneled fm index
    tfm_index<> tfm;
    load_from_file(tfm, argv[1]);

    auto *S = new uint32_t[tfm.size()];
    S[tfm.size() - 1] = 0;
    auto p = tfm.end();
    for (size_type i = 1; i < tfm.size(); i++) {
        S[tfm.size() - i - 1] = (uint32_t) tfm.backwardstep(p);
    }

    ofstream wf(argv[2], ios::out | ios::binary);
    for (int i = 0; i < tfm.size() - 1; i++) {
        wf.write((char *) &S[i], sizeof(uint32_t));
    }
}
