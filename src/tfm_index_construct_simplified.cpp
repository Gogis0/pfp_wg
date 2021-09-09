#include "functions.cpp"

int main() {
    // vector<string> dict = read_dict("../data/yeast.fasta.dict");
    // vector<uint> parse = read_parse("../data/yeast.fasta.parse");
    vector<uint> parse = {1, 2, 1, 2};  // aba aca aba aca
    vector<string> dict = {"a", "ab", "ac"};

    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    cout << wg.ordering << endl;
    cout << wg.dot_repr_ordered() << endl;
    wg_unparse(wg, dict);
    cout << wg.dot_repr_ordered() << endl;
    cout << (wg.is_valid()? "is valid" : "is invalid") << endl;


    return 0;
}

