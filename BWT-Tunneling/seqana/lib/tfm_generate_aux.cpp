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

using namespace std;
using namespace sdsl;

typedef typename sdsl::int_vector<>::size_type size_type;

void printUsage( char **argv ) {
	cerr << "USAGE: " << argv[0] << " TFMFILE OCCFILE" << endl;
	cerr << "TFMFILE:" << endl;
	cerr << "  File where to store the serialized trie" << endl;
};

int main(int argc, char **argv) {
	//check parameters
	if (argc != 3) {
		printUsage( argv );
		cerr << "At least 2 parameters expected" << endl;
		return 1;
	}
	
	//load tunneled fm index
	tfm_index<> tfm;
	load_from_file( tfm, argv[1] );

        uint32_t dphrases = 0;
	for (size_type i = 0; i < tfm.L.size(); i++) {
            uint32_t x = (uint32_t)tfm.L[i];
            if (x > dphrases) dphrases = x;
        }

        vector<vector<uint32_t>> phrase_sources(tfm.L.size());
        vector<uint32_t> phrase_occs(dphrases, (uint32_t)0);

	auto p = tfm.end();
        int last = 1;
	for (size_type i = 0; i < tfm.L.size(); i++) {
            //if (tfm.dout[i] == 1) last = i;
            auto c = (uint)tfm.L[i];
            //phrase_sources[last-1].push_back((uint)tfm.L[i]);
            phrase_occs[c]++;
            //cout << (uint)tfm.L[i] << " ";
	}
        // read the dictionary
        FILE *fdict = fopen(argv[2], "r");
        for (uint32_t i )

        // write the numbers of occurrences for each phrase
        FILE *focc = fopen(argv[3], "wb");
        for (uint32_t i = 0; i < dphrases; i++) {
            size_t s = fwrite(&phrase_occs[i], sizeof(phrase_occs[i]), 1, focc);
            if (s != 1) {
                cout << "Error writing to OCC file" << endl;
                exit(1);
            }
        }
        cout << "OCC file successfully written!" << endl;

        // write 
        /*
	for (size_type i=0; i<tfm.size()-1; i++) {
	    cout << S[i] << ' ';
	}
        */
}