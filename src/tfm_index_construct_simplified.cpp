#include "functions.cpp"

int main() {
    // vector<string> dict = read_dict("../data/yeast.fasta.dict");
    // vector<uint> parse = read_parse("../data/yeast.fasta.parse");
    vector<uint> parse = {1, 2, 1, 2};  // aba aca aba aca
    vector<string> dict = {"a", "ab", "ac"};

    tfm_index<> tfm;
    my_construct(tfm, parse);

    wheeler_graph wg = wheeler_graph(tfm);
    // cout << wg.dot_repr() << endl;
    wg_unparse(wg, dict);
    // cout << wg.dot_repr() << endl;

    wg.ordering = {0, 0, 0, 0, 0, 0, 0, 0};
    cout << wg.dot_repr() << endl;
    if (wg.is_valid())
        cout << "Is wheeler graph\n" << endl;
    else
        cout << "Is not wheeler graph\n" << endl;

    wg.ordering = {2, 5, 6, 7, 0, 3, 1, 4};

    cout << wg.dot_repr_ordered() << endl;
    if (wg.is_valid())
        cout << "Is wheeler graph\n" << endl;
    else
        cout << "Is not wheeler graph\n" << endl;

    return 0;
}

