#include "functions.cpp"

int main() {
    wheeler_graph wg = wheeler_graph(2);
    wg.add_edge(0, 1, 0, 0);
    wg.add_edge(0, 1, 0, 1);
    cout << wg.dot_repr() << endl;
    vector<string> dict = {"BABAABA"};
    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;

//    string T = "abcadeafgabcadeafga";
//    vector<string> E = {"a"};
//
//    vector<string> dict;
//    vector<uint> parse;
//     fill_dict_and_parse(T, E, dict, parse);
//
//    tfm_index<> tfm = tfm_create(parse);
//    cout << tfm_repr(tfm) << endl;
//    wheeler_graph wg = wheeler_graph(tfm);
//    wg_unparse(wg, dict);
//
//    cout << wg.dot_repr() << endl;

//    tfm_index<> tfm2 = wg_to_tfm(wg);
//    cout << tfm2.size() << endl;
//    cout << tfm_repr(tfm2) << endl;
//    for (unsigned char i : tfm2.L) {cout << to_string(i) << " ";}
//    cout << endl;
//    cout << tfm2.dout << endl;
//    cout << tfm2.din << endl;

    return 0;
}

