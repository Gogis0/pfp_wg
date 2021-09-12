#include "functions.cpp"

string reverse(const string &s) {
    string res(s.length(), 'x');
    for (uint i = 0; i < s.length(); i++) {
        res[s.length() - 1 - i] = s[i];
    }
    return res;
}

bool cmp(const pair<string, uint> &p1, const pair<string, uint> &p2) {
    string s1 = reverse(p1.first);
    string s2 = reverse(p2.first);
    return (s1.compare(s2));
}

int main() {
//    vector<string> dict = read_dict("../data/yeast.fasta.dict");
//    vector<uint> parse = read_parse("../data/yeast.fasta.parse");
//
//    dict.emplace_back("\u0002\u0002\u0002\u0002\u0002CCA");
//    vector<pair<string, uint>> v;
//    for (uint i = 0; i<dict.size(); i++) {
//        dict[i] = dict[i].substr(0, dict[i].length() - 4);
//        v.emplace_back(dict[i], i);
//    }
//    sort(v.begin(), v.end(), cmp);
//    vector<uint> v2;
//    for (uint i = 0; i < dict.size(); i++) {
//        dict[i] = v[i].first;
//        v2.emplace_back(v[i].second);
//    }
//    vector<uint> v3 = inverse_permutation(v2);
//    vector<uint> p = vector<uint> (parse.size(), 0);
//    for (uint i = 0; i < parse.size(); i++) {
//        p[i] = parse[v3[i]];
//    }
//    vector<string> dict = read_dict("../data/yeast.fasta.dict");
//    for (uint i = 0; i < dict.size(); i++) {
//        cout << dict[i] << endl;
//    }
    // T=$CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA$$$$
    vector<uint> parse = read_parse("../data/yeast.fasta.parse");
    vector<string> dict = {
                    "$$$$", //  7   -> 0
          "CACCCACACACACA", //  4   -> 1
            "CACCCACACACA", //  5   -> 2
                 "CACCACA", //  2   -> 3
                "CACCCACA", //  6   -> 4
                    "$CCA", //  1   -> 5
                   "CACCA", //  3   -> 6
    };

    cout << parse << endl;
    vector<uint> remap = {0, 5, 3, 6, 1, 2, 4, 0};
    for (uint i = 0; i < parse.size(); i++) {
        parse[i] = remap[parse[i]];
    }
    cout << parse << endl;

    tfm_index<> tfm;
    my_construct(tfm, parse);
    cout << dot_repr_tfm(tfm) << endl;
    wheeler_graph wg = wheeler_graph(tfm);
    cout << wg.dot_repr() << endl;
    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;
    cout << reverse(wg_string(wg)) << endl;
    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

