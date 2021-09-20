#include "functions.cpp"

int main() {
    vector<uint> parse = {1, 2, 3, 1, 2, 3};
    tfm_index<> tfm = tfm_create(parse);
    cout << tfm_repr(tfm) << endl;
    for (unsigned char i : tfm.L) {cout << to_string(i) << " ";}
    cout << endl;
    cout << tfm.dout << endl;
    cout << tfm.din << endl;
    cout << tfm.size() << endl;

    wheeler_graph wg = wheeler_graph(tfm);
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

