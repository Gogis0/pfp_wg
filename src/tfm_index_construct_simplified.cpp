#include "functions.cpp"

int main() {
    string text =
        "$CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACCACACC"
        "ACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA$$$$";
    vector<string> triggers = {"$CCA", "CACC", "$$$$"};

    vector<uint> parse;
    vector<string> dict;
    fill_dict_and_parse(text, triggers, dict, parse);
    cout << parse << endl;
    for (auto &str: dict) cout << str << endl;

    tfm_index<> tfm = tfm_create(parse);
    cout << tfm_repr(tfm) << endl;

    wheeler_graph wg = wheeler_graph(tfm);
    cout << wg.dot_repr() << endl;

    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;
    cout << wg_string(wg) << endl;

    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

