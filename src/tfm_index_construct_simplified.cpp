#include "functions.cpp"

string loadfile(const string &filename) {
    stringstream ss;
    ifstream file(filename);
    stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

int main() {
    string path = "../BWT-Tunneling/testdata/pizzachili/tmp";
    string T = "\x04\x04" + loadfile(path) + "\x03\x03";
    // cout << T << endl;

    vector<string> E = {"\x04\x04", "vf", "gc", "fk", "\x03\x03"};

    vector<string> dict;
    vector<uint> parse;
    fill_dict_and_parse(T, E, dict, parse);
    for (uint i = 0; i < dict.size(); i++) cout << dict[i] << endl;
    cout << parse << endl;

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

