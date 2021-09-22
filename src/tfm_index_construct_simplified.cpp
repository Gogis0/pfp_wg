#include "functions.cpp"

int main() {
//    string T = "abcadeafgahiajkabcadeafgahiajka";
    string T = "abacadabacada";
    vector<string> E = {"a"};

    vector<string> dict;
    vector<uint> parse;
    fill_dict_and_parse(T, E, dict, parse);
    cout << parse << endl;
    for (uint i = 0; i < dict.size(); i++) cout << dict[i] << endl;

    tfm_index<> tfm = tfm_create(parse);
    cout << "TFM1\n" << tfm_repr(tfm) << endl;
    cout << tfm.size() << endl;
    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    position_end(wg, dict[0]);
    wg_find_ordering(wg);

    cout << wg.dot_repr() << endl;

    tfm_index<> tfm2 = wg_to_tfm(wg);
    cout << tfm2.size() << endl;
    cout << "TFM2\n" << tfm_repr(tfm2) << endl;
    for (unsigned char i : tfm2.L) {cout << to_string(i) << " ";}
    cout << endl;
    cout << tfm2.dout << endl;
    cout << tfm2.din << endl;

    tfm_index<> tfm3;
    wt_blcd<> wt = construct_from_il({4, 4, 1, 2, 3, 1, 1, 1, 1});
    bit_vector dout = int_vector<1>({1, 1, 1, 1, 1, 1, 0, 1, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 1, 1, 1, 1, 0, 1});
    uint64_t text_len = 13;

    construct_tfm_index_tmp(tfm3, text_len, move(wt), move(dout), move(din));
    cout << "TFM3\n" << tfm_repr(tfm3) << endl;

    print_order(wg);

    return 0;
}

