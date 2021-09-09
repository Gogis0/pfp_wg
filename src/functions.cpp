#include <iostream>
#include <utility>
#include <sdsl/util.hpp>
#include "tfm_index.hpp"

using namespace std;
using namespace sdsl;

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

    string dot_repr() const {
        stringstream ss;
        ss  << "\t" << u
            << " -> " << v
            << " [label = \"t=" << t
            << "\\nl=" << l
            << "\"]\n";
        return ss.str();
    }
};

uint my_select(const vector<uint> &v, uint elem, uint tier) {
    for (uint i = 0; i < v.size(); i++) {
        if (v[i] == elem) {
            if (tier == 0) return i;
            else tier--;
        }
    }
    return UINT32_MAX;
}

uint my_rank(const vector<uint> &v, uint elem, uint pos) {
    uint rank = 0;
    for (uint i = 0; i < pos; i++) {
        if (v[i] == elem) rank++;
    }
    return rank;
}

class wheeler_graph{
public:
    uint n_vertices;
    uint n_edges;
    vector<uint> source_nodes;
    vector<uint> destination_nodes;
    vector<uint> ordering;
    vector<uint> labels;
    vector<uint> tunnel_num;

    explicit wheeler_graph(uint V) {
        n_vertices = V;
        n_edges = 0;
    }

    explicit wheeler_graph(const tfm_index<> &tfm) {
        n_vertices = tfm.L.size();
        for (uint i = 0; i < n_vertices; i++) {ordering.push_back(i);}

        n_edges = tfm.din.size();
        edge e = edge(0, 0, 0, 0);
        for (auto i = tfm.size(); i > 0; i--) {
            e = e.get_next(tfm);
            source_nodes.push_back(e.u);
            destination_nodes.push_back(e.v);
            labels.push_back(e.l);
            tunnel_num.push_back(e.t);
        }

    }

    pair<uint, uint> end(){
        return {0, 0};
    };

    void backward(pair<uint, uint> &pos) {
        uint i = my_select(source_nodes, pos.first, pos.second); // from (source_nodes) select (pos.second)th (pos.first)
        pos.first = destination_nodes[i];
        pos.second = tunnel_num[i];
    };

    void forward(pair<uint, uint> &pos) {
        uint i = my_select(destination_nodes, pos.first, pos.second);
        pos.first = source_nodes[i];
        pos.second = my_rank(source_nodes, pos.first, i);
    };

    uint preceding_letter(const pair<uint,uint> &pos) {
        uint i = my_select(destination_nodes, pos.first, pos.second);
        return labels[i];
    }

    bool is_valid() {
        for (uint i = 0; i < n_edges; i++) {
            for (uint j = 0; j < n_edges; j++) {
                if (labels[i] < labels[j]
                    && ordering[destination_nodes[i]] >= ordering[destination_nodes[j]]
                ) {
                    return false;
                }
                if (labels[i] == labels[j]
                    && ordering[source_nodes[i]] < ordering[source_nodes[j]]
                    && ordering[destination_nodes[i]] > ordering[destination_nodes[j]]
                ) {
                    return false;
                }
            }
        }
        return true;
    }

    void add_edge(uint source, uint destination, uint label, uint tunnel) {
        source_nodes.push_back(source);
        destination_nodes.push_back(destination);
        labels.push_back(label);
        tunnel_num.push_back(tunnel);
        n_edges++;
    }

    void remove_edge(uint edge) {
        source_nodes.erase(source_nodes.begin() + edge);
        destination_nodes.erase(destination_nodes.begin() + edge);
        labels.erase(labels.begin() + edge);
        tunnel_num.erase(tunnel_num.begin() + edge);
        n_edges--;
    }

    string dot_repr() {
        stringstream ss;
        ss << "digraph G {\n";
        for (uint i = 0; i < n_edges; i++) {
            ss  << "\t" << source_nodes[i]
                << " -> " << destination_nodes[i]
                << " [label = \"t=" << tunnel_num[i]
                << "\\nl=" << labels[i]
                << "\"]\n";
        }
        ss << "}";
        return ss.str();
    }

    string dot_repr_ordered() {
        stringstream ss;
        ss << "digraph G {\n";
        for (uint i = 0; i < n_edges; i++) {
            ss  << "\t" << ordering[source_nodes[i]]
                << " -> " << ordering[destination_nodes[i]]
                << " [label = \"t=" << tunnel_num[i]
                << "\\nl=" << labels[i]
                << "\"]\n";
        }
        ss << "}";
        return ss.str();
    }
};

void expand_edge(wheeler_graph &wg, uint edge, const string &label) {
    uint end = wg.destination_nodes[edge];
    uint current = wg.source_nodes[edge];
    uint i;
    for (i = 0; i < label.length()-1; i++) {
        uint new_vertex = wg.n_vertices;
        wg.n_vertices++;
        wg.add_edge(current, new_vertex, label[i], wg.tunnel_num[i]);
        current = new_vertex;
    }
    wg.add_edge(current, end, label[i], wg.tunnel_num[i]);
    wg.remove_edge(edge);
}

int cmp_vertices(wheeler_graph &wg, pair<uint, uint> v1, pair<uint, uint> v2, vector<uint> original_ordering) {
    // precondition: vertices are not equal
    // v1 > v2 => 1, v1 < v2 => -1
    if (v1.first < original_ordering.size() && v2.first < original_ordering.size()) {
        if (original_ordering[v1.first] < original_ordering[v2.first]) {
            return -1;
        } else if (original_ordering[v2.first] < original_ordering[v1.first]) {
            return 1;
        }
    }

    if (wg.preceding_letter(v1) < wg.preceding_letter(v2)) {
        return -1;
    } else if (wg.preceding_letter(v2) < wg.preceding_letter(v1)) {
        return 1;
    } else {
        wg.forward(v1);
        wg.forward(v2);
        return cmp_vertices(wg, v1, v2, original_ordering);
    }
}

vector<uint>::iterator partition(
    vector<uint>::iterator start, vector<uint>::iterator end, wheeler_graph &wg, const vector<uint> &orig_ord
) {
    auto p = (start + (end - start)/2);
    auto pivot = pair<uint, uint> {*p, 0};
    auto i = start;

    for (auto j = start; j < end; j++) {
        if (j == p) continue;
        auto cur_vert = pair<uint, uint> {*j, 0};
        // if cur_vert <= pivot:
        if (cmp_vertices(wg, cur_vert, pivot, orig_ord) == -1) {
            swap(*i, *j);
            i++;
        }
    }
    swap(*p, *i);
    return i;
}

void my_sort(
    vector<uint>::iterator start, vector<uint>::iterator end, wheeler_graph &wg, const vector<uint> &orig_ord
) {
    if ((end - start) <= 1) return;
    auto middle = partition(start, end, wg, orig_ord);
    my_sort(start, middle, wg, orig_ord);
    my_sort(middle + 1, end, wg, orig_ord);
}

vector<uint> inverse_permutation(const vector<uint> &v) {
    vector<uint> res(v.size(), 0);
    for (uint i = 0; i < v.size(); i++) {
        res[v[i]] = i;
    }
    return res;
}

void wg_find_ordering(wheeler_graph &wg) {
    vector<uint> original_ordering = wg.ordering;
    vector<uint> new_ordering;
    for (uint i = 0; i < wg.n_vertices; i++) {new_ordering.push_back(i);}
    my_sort(new_ordering.begin(), new_ordering.end(), wg, original_ordering);
    wg.ordering = inverse_permutation(new_ordering);
}

void wg_unparse(wheeler_graph &wg, vector<string> &dict) {
    uint n_edges = wg.n_edges;
    for (uint i = 0; i < n_edges; i++) {
        expand_edge(wg, 0, dict[wg.labels[0]]);
    }
    wg_find_ordering(wg);
}

vector<uint> read_parse(const string &infile) {
    vector<uint> v;
    FILE *parse = fopen(infile.c_str(), "r");

    fseek(parse, 0, SEEK_END);
    size_t n;
    n = ftell(parse) / sizeof(int);
    fseek(parse, 0, SEEK_SET);

    int tmp;
    for (int i = 0; i < n; i++) {
        fread(&tmp, sizeof(int), 1, parse);
        v.push_back(tmp);
    }
    return v;
}

vector<string> read_dict(const string &infile) {
    vector<string> v;
    FILE *dict = fopen(infile.c_str(), "rb");

    fseek(dict, 0, SEEK_END);
    size_t n = ftell(dict) / sizeof (char);
    fseek(dict, 0, SEEK_SET);

    char c;
    string str;
    for (uint i = 0; i < n; i++) {
        fread(&c, sizeof (char), 1, dict);
        if (c == 0x01){
            v.push_back(str);
            str = "";
        } else {
            str += c;
        }
    }
    return v;
}

int_vector<> init_parse(const vector<uint> &il) {
    auto res = int_vector<>(il.size(), 1, 32);
    uint i = 0;
    for (auto value : il) {
        res[i] = value;
        i++;
    }
    return res;
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

string dot_repr_tfm(const tfm_index<> &tfm) {
    edge e = edge(0, 0, 0, 0);
    stringstream ss;
    ss << "digraph G {\n";
    for (auto i = tfm.size(); i > 0; i--) {
        e = e.get_next(tfm);
        ss << e.dot_repr();
    }
    ss << "}";
    return ss.str();
}

void my_construct(tfm_index<> &tfm,  vector<uint> &parse_vec) {
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    cache_config config = cache_config(false, "../data/", "tmp");

    int_vector<> parse = init_parse(parse_vec);
    csa_wt<wt_blcd_int<>> csa;
    construct_im(csa, parse);

    //find minimal edge-reduced DBG and store kmer bounds in a bitvector B
    bit_vector B;
    dbg_algorithms::find_min_dbg(csa, B, config);
    // cout << "B: " << B << endl;

    //use bitvector to determine prefix intervals to be tunneled
    bit_vector dout = B;
    bit_vector din = B;
    dbg_algorithms::mark_prefix_intervals(csa, dout, din);
    // cout << "din:\t" << din << "\ndout:\t" << dout << "\n\n";

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
