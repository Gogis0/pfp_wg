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
    cout << "I:\t" << tfm.din << endl;              // I[1..|v|]
    cout << "O:\t" << tfm.dout << endl;             // O[1..|v|]
    cout << "C:\t";
    for (uint i = 0; i <= tfm.L.sigma; i++) {
        cout << tfm.C[i] <<  " " ;                  // C[1..|A|]
    }
    cout << endl << endl;
}

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

void my_construct(tfm_index<> &tfm, int_vector<> &parse) {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);

    //find minimal edge-reduced DBG and store kmer bounds in a bitvector B
    bit_vector B;
    dbg_algorithms::find_min_dbg(csa, B, config);
    cout << "B: " << B << endl;

    //use bitvector to determine prefix intervals to be tunneled
    bit_vector dout = B;
    bit_vector din = B;
    dbg_algorithms::mark_prefix_intervals(csa, dout, din);
    cout << "din:\t" << din << "\ndout:\t" << dout << "\n\n";

    //create a buffer for newly constructed L
    string tmp_file_name = "tmp_123";
    int_vector_buffer<8> L_buf(tmp_file_name, std::ios::out);

    for (ulong i = 0; i < csa.size(); i++) {
        if (din[i] == 1) L_buf.push_back(csa.wavelet_tree[i]);
    }

    //remove redundant entries from L, dout and din
    ulong p = 0;
    ulong q = 0;
    for (ulong i = 0; i < csa.size(); i++) {
        // cout << "din:\t" << din << "\ndout:\t" << dout << "\n";
        if (din[i] == 1) {
            dout[p] = dout[i];
            // cout << "dout " << p << "<-" << i << endl;
            p++;
        }
        if (dout[i] == 1) {
            din[q] = din[i];
            // cout << "din " << q << "<-" << i << endl;
            q++;
        }
    }
    // cout << "din:\t" << din << "\ndout:\t" << dout << "\n";
    dout[p] = 1;
    p++;
    din[q] = 1;
    q++;
    dout.resize(p);
    din.resize(q);
    // cout << "din:\t" << din << "\ndout:\t" << dout << "\n";

    uint64_t text_len = csa.size();

    construct_tfm_index(tfm, text_len, std::move(L_buf), std::move(dout), std::move(din));
    //remove buffer for L
    sdsl::remove(tmp_file_name);
    // cout << "\n" << endl;
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

