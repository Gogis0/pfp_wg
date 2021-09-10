#include "functions.cpp"

int main() {
    // T=acbdacbda, E={a, b}
    // acb bda acb bda aa -> ac bd a -> a ac bd
    vector<string> dict = {"a", "ac", "bd"};
    vector<uint> parse = {1, 2, 1, 2};  // 0 at the end is implied

    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;
    cout << wg_string(wg) << endl;
    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

