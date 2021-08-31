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
    cout << "L:\t";
    for (unsigned char i : tfm.L) {
        cout << to_string(i) << " ";                // L[1..|e|]
    }
    cout << endl;
    cout << "O:\t" << tfm.dout << endl;             // O[1..|v|]
    cout << "I:\t" << tfm.din << endl;              // I[1..|v|]
    cout << "C:\t";
    for (uint i = 0; i <= tfm.L.sigma; i++) {
        cout << tfm.C[i] <<  " " ;                  // C[1..|A|]
    }
    cout << endl << endl;
}

// next step: read BWT-Tunneling/seqana/include/dbg_algorithms.hpp

class edge{
public:
    ulong u;    // in node
    ulong v;    // out node
    ulong l;    // label
    ulong t;    // tunnel number

    edge() = default;
    edge(uint in_node, uint out_node, uint label, uint tunnel){
        u = in_node;
        v = out_node;
        l = label;
        t = tunnel;
    }
    edge get_next(const tfm_index<> &tfm) const {
        auto L = tfm.L;
        auto O = tfm.dout;
        auto I = tfm.din;
        auto C = tfm.C;
        auto rI = tfm.din_rank;
        auto sI = tfm.din_select;
        auto sO = tfm.dout_select;

        edge next = edge();
        next.u = v;
        next.t = t;
        next.l = L.inverse_select(next.u).second;
        next.v = L.inverse_select(next.u).first + C[next.l];
        if (I[next.v] == 0) {
            next.t = next.v - sI(rI(next.v + 1));
        }
        next.v = sO(rI(next.v + 1));
        if (O[next.v + 1] == 0) {
            next.v += next.t; next.t = 0;
        }
        return next;
    }
    void print() const {
        cout << "Tunnel " << t << " from " << u << " to " << v << " with label " << l << ".\n";
    }
};


void print_original(const tfm_index<> &tfm) {
    edge e = edge(0, 0, 0, 0);
    for (auto i = tfm.size(); i > 0; i--) {
        e = e.get_next(tfm);
        e.print();
    }
    cout << endl;
}

wt_blcd<> construct_from_vector(initializer_list<uint32_t> il);

int main() {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    tfm_index<> tfm;
    auto parse = init_parse({1, 2, 1, 2});
    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);
    construct_tfm_index(tfm, move(csa), config);
    print_tfm(tfm);
    print_original(tfm);



    tfm_index<> tfm2;
    wt_blcd<> wt = construct_from_vector({2, 2, 0, 1});
    bit_vector dout = int_vector<1>({1, 1, 0, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 0, 1});
    uint64_t text_len = tfm.size();

    construct_tfm_index_tmp(tfm2, text_len, move(wt), move(dout), move(din));
    print_tfm(tfm2);
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
