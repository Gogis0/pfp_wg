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
            << "\t-> " << v
            << "\t[label = \"t=" << t
            << "\\nl=" << l
            << "\"]\n";
        return ss.str();
    }
};

string reverse(const string &s) {
    string res(s.length(), 'x');
    for (uint i = 0; i < s.length(); i++) {
        res[s.length() - 1 - i] = s[i];
    }
    return res;
}

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
    pair<uint, uint> end = {0, 0};

    explicit wheeler_graph(uint V) {
        n_vertices = V;
        n_edges = 0;
    }

    explicit wheeler_graph(const tfm_index<> &tfm) {
        n_vertices = tfm.L.size();
        for (uint i = 0; i < n_vertices; i++) {ordering.push_back(i);}

        // n_edges = tfm.din.size();
        n_edges = tfm.size();
        edge e = edge(0, 0, 0, 0);
        for (auto i = tfm.size(); i > 0; i--) {
            e = e.get_next(tfm);
            source_nodes.push_back(e.u);
            destination_nodes.push_back(e.v);
            labels.push_back(e.l);
            tunnel_num.push_back(e.t);
        }

    }

    pair<uint, uint> get_end() const {
        return end;
    };

    void backward(pair<uint, uint> &pos) const {
        uint i = my_select(source_nodes, pos.first, pos.second); // from (source_nodes) select (pos.second)th (pos.first)
        pos.first = destination_nodes[i];
        pos.second = tunnel_num[i];
    };

    void forward(pair<uint, uint> &pos) const {
        uint i = my_select(destination_nodes, pos.first, pos.second);
        pos.first = source_nodes[i];
        pos.second = my_rank(source_nodes, pos.first, i);
    };

    uint preceding_letter(const pair<uint,uint> &pos) const {
        uint i = my_select(destination_nodes, pos.first, 0);
        return labels[i];
        // return labels[labels.size() - 1 - i];
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
            ss  << "\t\"" << source_nodes[i];
            if (source_nodes[i] < ordering.size()) ss << ":" << ordering[source_nodes[i]];
            ss  << "\"\t-> \"" << destination_nodes[i];
            if (source_nodes[i] < ordering.size()) ss << ":" << ordering[destination_nodes[i]];
            ss  << "\"\t[label = \"t=" << tunnel_num[i]
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
        wg.add_edge(current, new_vertex, label[label.length() - 1 - i], wg.tunnel_num[edge]);
        current = new_vertex;
    }
    wg.add_edge(current, end, label[0], wg.tunnel_num[edge]);  // (label.length() - 1 - i) == 0 at this point
}

void expand_tunneled_edge(wheeler_graph &wg, uint edge, const string &label) {
    uint source = wg.source_nodes[edge];
    uint destination = wg.destination_nodes[edge];
    uint tun = wg.tunnel_num[edge];
    wg.remove_edge(edge);

    pair<uint, uint> p = {destination, 0};
    uint current = p.first;
    for (uint i = 0; i < label.length() - 1; i++) {
        wg.forward(p);
        wg.add_edge(p.first, current, label[i], tun);
        current = p.first;
    }
    wg.add_edge(source, current, label[label.length() - 1], tun);
}

int cmp_vertices(wheeler_graph &wg, pair<uint, uint> v1, pair<uint, uint> v2, vector<uint> original_ordering) {
    // precondition: vertices are not equal
    uint i, j;
    i = wg.preceding_letter(v1);
    j = wg.preceding_letter(v2);
    if (i < j) { return -1; }
    if (i > j) { return 1; }

    uint n = original_ordering.size();
    if (v1.first < n && v2.first < n) {
        i = original_ordering[v1.first];
        j = original_ordering[v2.first];
        if (i < j) { return -1; }
        if (i > j) { return 1; }
    }

    wg.forward(v1);
    wg.forward(v2);
    return cmp_vertices(wg, v1, v2, original_ordering);
}

vector<uint>::iterator partition(
    vector<uint>::iterator start, vector<uint>::iterator end, wheeler_graph &wg, const vector<uint> &orig_ord
) {
    auto p = end - 1;
    auto pivot = pair<uint, uint> {*p, 0};
    auto i = start;

    for (auto j = start; j < p; j++) {
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
    vector<uint>::iterator start,   // first element
    vector<uint>::iterator end,     // past_the_end element
    wheeler_graph &wg,
    const vector<uint> &orig_ord
) {
    // for (auto it = start; it < end; it++) cout << *it << "\t";
    // cout << endl;
    if ((end - start) <= 1) return;
    auto middle = partition(start, end, wg, orig_ord);
    // cout << "middle:" << *middle << endl;
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
    // for (uint i = 0; i < wg.n_vertices; i++) {wg.ordering[i] = wg.n_vertices - wg.ordering[i];}
}

void position_end(wheeler_graph &wg, const string &end_label) {
    auto end = wg.get_end();
    for (uint i = 0; i < end_label.length()-1; i++)
        wg.forward(end);
    wg.end = end;
}

void wg_unparse(wheeler_graph &wg, vector<string> &dict) {
    uint n_edges = wg.n_edges;
    string label;
    for (uint i = 0; i < n_edges; i++) {
        label = dict[wg.labels[0]];
        if (wg.tunnel_num[0] == 0) {
            expand_edge(wg, 0, label);
            wg.remove_edge(0);
        } else {
            expand_tunneled_edge(wg, 0, label);
        }
    }
}

string wg_string(const wheeler_graph &wg) {
    stringstream ss;
    auto pos = wg.get_end();
    char c = (char)wg.preceding_letter(pos);
    ss << c;
    for (uint i = 0; i < wg.n_edges - 1; i++) {
        wg.backward(pos);
        c = (char)wg.preceding_letter(pos);
        ss << c;
    }
    return reverse(ss.str());
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

wt_blcd<> construct_from_il(initializer_list<uint32_t> il) {
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

wt_blcd<> construct_from_vector(const vector<uint> &v) {
    wt_blcd<> wt;
    const uint8_t size = wt_blcd<>::alphabet_category::WIDTH;
    int_vector<size> text;
    text.resize(v.size());
    for (uint i = 0; i < v.size(); i++) {
        text[i] = v[i];
    }

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

string tfm_repr(const tfm_index<> &tfm) {
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

tfm_index<> tfm_create(vector<uint> &parse_vec) {
    tfm_index<> tfm;
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
    return tfm;
}

bool cmp(const pair<string, uint> &p1, const pair<string, uint> &p2) {
    string s1 = p1.first;
    string s2 = p2.first;
    return s1 < s2;
}

void create_parse(const string &text, const vector<string> &triggers, map<string, uint> &dict, vector<uint> &parse) {
    uint w = triggers[0].length();
    uint p = 0;
    uint phrase_start = 0;
    for (uint i = phrase_start + 1; i < text.length(); i++) {
        string substr = text.substr(i, w);
        for (const auto & trigger : triggers) {
            if (substr == trigger) {
                string phrase = text.substr(phrase_start, i + w - phrase_start);
                phrase_start = i;
                auto it = dict.find(phrase);
                if (it == dict.end()) {
                    dict.insert({phrase, p});
                    parse.push_back(p);
                    p++;
                } else {
                    parse.push_back(it->second);
                }
            }
        }
    }
    dict.insert({triggers[triggers.size() - 1] + triggers[0], p});
    parse.push_back(p);
}

void fill_dict_and_parse(const string &text, const vector<string> &E, vector<string> &dict, vector<uint> &parse){
    map<string, uint> D;
    vector<uint> P;

    create_parse(text, E, D, P);

    // sort dict
    vector<pair<string, uint>> v;
    v.reserve(D.size());
    for (const auto & pair : D) {
        v.emplace_back(pair.first, pair.second);
    }
    sort(v.begin(), v.end(), cmp);

    vector<uint> remap = vector<uint>(P.size(), 0);
    for (uint i = 0; i < D.size(); i++) {
        dict.push_back(v[i].first);
        remap[v[i].second] = i;
    }

    parse = vector<uint>(P.size(), 0);
    for (uint i = 0; i < P.size(); i++) {
        parse[i] = remap[P[i]];
    }

    uint w = E[0].length();
    for (auto & word : dict) { word = word.substr(0, word.length() - w); }
    parse = vector<uint>(parse.begin(), parse.end() - 1);
}

uint count_out_nodes(const wheeler_graph &wg, uint node) {
    set<uint> s;
    for (uint i=0; i < wg.n_edges; i++) {
        if (wg.source_nodes[i] == node) {
            s.insert(wg.destination_nodes[i]);
        }
    }
    return s.size();
}

uint count_in_nodes(const wheeler_graph &wg, uint node) {
    set<uint> s;
    for (uint i=0; i < wg.n_edges; i++) {
        if (wg.destination_nodes[i] == node) {
            s.insert(wg.source_nodes[i]);
        }
    }
    return s.size();
}

tfm_index<> wg_to_tfm(const wheeler_graph &wg) {
    vector<uint> labels;

    vector<uint> inv_order = inverse_permutation(wg.ordering);
    for (uint i = 0; i < inv_order.size(); i++) {
        pair<uint, uint> p(inv_order[i], 0);
        wg.backward(p);
        labels.push_back(wg.preceding_letter(p));
    }

    bit_vector din(wg.n_edges, 0);
    bit_vector dout(wg.n_edges, 0);

    uint i = 0;
    uint j = 0;
    for (uint p = 0; p < wg.n_vertices; p++) {
        dout[i] = 1;
        din[j] = 1;

        pair<uint, uint> pos(p, 0);
        wg.forward(pos);
        uint p_out_nodes = count_out_nodes(wg, pos.first);
        i += p_out_nodes;

        uint p_in_nodes = count_in_nodes(wg, p);
        j += p_in_nodes;

        p += p_out_nodes - 1;
    }
    dout[i] = 1;
    din[j] = 1;
    i++;
    j++;
    dout.resize(i);
    din.resize(j);

    tfm_index<> tfm;
    wt_blcd<> wt1 = construct_from_vector(labels);
    // wt_blcd<> wt2 = construct_from_il({3, 3, 0, 1, 2});
    construct_tfm_index_tmp(tfm, wg.n_edges, move(wt1), move(dout), move(din));
    return tfm;
}
