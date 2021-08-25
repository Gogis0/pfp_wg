#include <iostream>
#include <utility>
#include <sdsl/util.hpp>
#include "../BWT-Tunneling/seqana/include/tfm_index.hpp"

using namespace std;
using namespace sdsl;

int_vector<> init_vector(initializer_list<uint32_t> il) {
    auto res = int_vector<>(il.size(), 1, 32);
    uint i = 0;
    for (auto value : il) {
        res[i] = value;
        i++;
    }
    return res;
}

int main() {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    tfm_index<> tfm;
    // auto parse = init_vector({1, 2, 3, 1, 2, 3, 4});
    auto parse = init_vector({1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4});
    cout << util::to_string(parse) << endl;
    cout << to_string(parse.width()) << endl;
    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);
    construct_tfm_index(tfm, move(csa), config);
    // construct_tfm_index(tfm, , L_buf, dout din);

    for (unsigned char i : tfm.L) {
        cout << to_string(i) << " ";    // L[1..|e|]
    }
    cout << endl;
    cout << tfm.dout << endl;           // O[1..|v|]
    cout << tfm.din << endl;            // I[1..|v|]
    cout << tfm.C << endl;              // C[1..|A|]

    store_to_file(tfm, "../data/yeast.tunnel");
    return 0;
}