#include "functions.cpp"

int main() {
    //              0  1  1  1   1  1  2  3  3
    vector<uint> v = {1, 2, 2, 3, 41, 1, 1, 2, 5};
    wt_blcd<> wt = construct_from_vector(v);
    wt_blcd_int<> wt2;
    for (uint i = 0; i < v.size(); i++)
        cout << wt.rank(i, 1) << " ";
    cout << endl;

    for (uint i = 1; i <= wt.rank(wt.size(), 1); i++)
        cout << wt.select(i, 1) << " ";
    cout << endl;

    return 0;
}
