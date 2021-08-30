#include <iostream>
#include <utility>
#include <sdsl/util.hpp>
#include "../BWT-Tunneling/seqana/include/tfm_index.hpp"

using namespace std;
using namespace sdsl;

int_vector<> init_parse(initializer_list<uint32_t> il) {
    auto res = int_vector<>(il.size(), 1, 32);
    uint i = 0;
    for (auto value : il) {
        res[i] = value;
        i++;
    }
    return res;
}

void print_tfm(const tfm_index<> &tfm) {
    for (unsigned char i : tfm.L) {
        cout << to_string(i) << " ";    // L[1..|e|]
    }
    cout << endl;
    cout << tfm.dout << endl;           // O[1..|v|]
    cout << tfm.din << endl;            // I[1..|v|]
    for (uint i = 0; i <= tfm.L.sigma; i++) {
        cout << tfm.C[i] <<  " " ;      // C[1..|A|]
    }
    cout << endl << endl;
}

wt_blcd<> construct_from_vector(initializer_list<uint32_t> il);

int main() {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    tfm_index<> tfm;
    // auto parse = init_parse({1, 2, 3, 1, 2, 3, 4});
    // auto parse = init_parse({1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4});
    auto parse = init_parse({1, 2, 1, 2});
    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);
    construct_tfm_index(tfm, move(csa), config);
    print_tfm(tfm);



    tfm_index<> tfm2;
    wt_blcd<> wt = construct_from_vector({4, 0, 3, 3, 3, 1, 2, 3});
    bit_vector dout = int_vector<1>({1, 1, 0, 0, 0, 1, 1, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 0, 0, 0, 1, 1});
    uint64_t text_len = tfm.size();

    construct_tfm_index_tmp(tfm2, text_len, move(wt), move(dout), move(din));
    // print_tfm(tfm2);

    store_to_file(tfm, "../data/yeast.tunnel");
    store_to_file(tfm2, "../data/yeast.tunnel2");
    return 0;
}

wt_blcd<> construct_from_vector(initializer_list<uint32_t> il) {
    wt_blcd<> wt;
    const uint8_t size = wt_blcd<>::alphabet_category::WIDTH;
    int_vector<size> text = int_vector<size>(il);

    string tmp_file_name = "tmp";
    store_to_file(text, tmp_file_name);

    int_vector_buffer<size> text_buf(tmp_file_name);
    wt_blcd<> tmp(text_buf, text_buf.size());
    wt.swap(tmp);
    remove(tmp_file_name);
    return wt;
}
