#include "functions.cpp"

int main() {
    // T=acbdacbda, E={a, b}
    string T = "acbdacbda";
    vector<string> E = {"b", "a"};

    vector<string> dict;
    vector<uint> full_parse;
    fill_dict_and_parse(T, E, dict, full_parse);
    cout << full_parse << endl;
    for (auto &str: dict) cout << str << endl;

    tfm_index<> tfm;
    vector<uint> parse = vector<uint>(full_parse.begin(), full_parse.end() - 1);
    cout << parse << endl;
    my_construct(tfm, parse);
    cout << dot_repr_tfm(tfm) << endl;
    wheeler_graph wg = wheeler_graph(tfm);
    cout << wg.dot_repr() << endl;
    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;
    cout << wg_string(wg) << endl;
    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

