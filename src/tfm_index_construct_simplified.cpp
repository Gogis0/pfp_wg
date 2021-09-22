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
    cout << tfm_repr(tfm) << endl;
    cout << tfm.size() << endl;
    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    position_end(wg, dict[0]);
    wg_find_ordering(wg);

    cout << wg.dot_repr() << endl;

    tfm_index<> tfm2 = wg_to_tfm(wg);
    cout << tfm2.size() << endl;
    cout << tfm_repr(tfm2) << endl;
    for (unsigned char i : tfm2.L) {cout << to_string(i) << " ";}
    cout << endl;
    cout << tfm2.dout << endl;
    cout << tfm2.din << endl;

    return 0;
}

