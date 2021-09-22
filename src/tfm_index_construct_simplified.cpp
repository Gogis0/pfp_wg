#include "functions.cpp"

int main() {
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

    return 0;
}

