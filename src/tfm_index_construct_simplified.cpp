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
    cout << tfm.L.size() <<  endl;
    for (uint i = 0; i < tfm.L.size(); i++) {   // How to get alphabet size?
        cout << tfm.C[i] <<  " " ;
    }
    cout << endl;
    // cout << tfm.C << endl;              // C[1..|A|]
}

wt_blcd<> construct_from_vector(initializer_list<uint32_t> il);

int main() {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    tfm_index<> tfm;
    // auto parse = init_parse({1, 2, 3, 1, 2, 3, 4});
    auto parse = init_parse({1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4});
    // cout << util::to_string(parse) << endl;
    // cout << to_string(parse.width()) << endl;
    // cout << parse << endl;
    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);
    construct_tfm_index(tfm, move(csa), config);
    print_tfm(tfm);

    tfm_index<> tfm2;
    // wt_blcd<> wt = tfm.L;
    wt_blcd<> wt = construct_from_vector({4, 0, 3, 3, 3, 1, 2, 3});
    // int_vector_buffer<>& buf = int_vector_buffer<>()
    // wt_blcd_int<> = wt_blcd_int();
    bit_vector dout = int_vector<1>({1, 1, 0, 0, 0, 1, 1, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 0, 0, 0, 1, 1});
    uint64_t text_len = tfm.size();

    construct_tfm_index_tmp(tfm2, text_len, move(wt), move(dout), move(din));
    print_tfm(tfm2);

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

wt_huff<> construct_wt() {
    string file = "@testinput.iv8";
    wt_huff<> wt;
    uint8_t num_bytes = 1;
    // construct(wt, file, num_bytes);
    // void construct(t_index& idx, std::string file, uint8_t num_bytes=0)
    tMSS file_map;
    cache_config config;
    if (is_ram_file(file)) {
        config.dir = "@";
    }
    // construct(wt, file, config, num_bytes);
    // void construct(t_index& idx, const std::string& file, cache_config& config, uint8_t num_bytes=0)
    typename wt_huff<>::index_category index_tag;
    // construct(wt, file, config, num_bytes, index_tag);
    // void construct(t_index& idx, const std::string& file, cache_config& config, uint8_t num_bytes, wt_tag)
    auto event = memory_monitor::event("construct wavelet tree");
    if ((wt_huff<>::alphabet_category::WIDTH==8 and num_bytes <= 1)
        or (wt_huff<>::alphabet_category::WIDTH==0 and num_bytes != 'd')) {
        int_vector_buffer<wt_huff<>::alphabet_category::WIDTH> text_buf(file, std::ios::in, 1024*1024, num_bytes*8, (bool)num_bytes);
        wt_huff<> tmp(text_buf, text_buf.size());
        wt.swap(tmp);
    } else {
        int_vector<wt_huff<>::alphabet_category::WIDTH> text;
        load_vector_from_file(text, file, num_bytes);
        std::string tmp_key = util::to_string(util::pid())+"_"+util::to_string(util::id());
        std::string tmp_file_name = cache_file_name(tmp_key, config);
        store_to_file(text, tmp_file_name);
        util::clear(text);
        {
            int_vector_buffer<wt_huff<>::alphabet_category::WIDTH> text_buf(tmp_file_name);
            wt_huff<> tmp(text_buf, text_buf.size());
            wt.swap(tmp);
        }
        sdsl::remove(tmp_file_name);
    }
    return wt;
}