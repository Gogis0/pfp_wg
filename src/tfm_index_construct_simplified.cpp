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
    return s1 < s2;
}

//    vector<string> dict = read_dict("../data/yeast.fasta.dict");
//    vector<uint> parse = read_parse("../data/yeast.fasta.parse");


void create_parse(const string &T, const vector<string> &E, map<string, uint> &D, vector<uint> &P) {
    uint p = 0;
    uint phrase_start = 0;
    for (uint i = phrase_start + 1; i < T.length(); i++) {
        string tmp = T.substr(i, 4);
        for (const auto & j : E) {
            if (tmp == j) {
                string phrase = T.substr(phrase_start, i - phrase_start);
                phrase_start = i;
                auto it = D.find(phrase);
                if (it == D.end()) {
                    D.insert({phrase, p});
                    P.push_back(p);
                    p++;
                } else {
                    P.push_back(it->second);
                }
            }
        }
    }
    D.insert({E[E.size()-1], p});
    P.push_back(p);
}

int main() {
    string T = "$CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA$$$$";
    vector<string> E = {"$CCA", "CACC", "$$$$"};

    map<string, uint> D;
    vector<uint> P;

    create_parse(T, E, D, P);

    // sort dict
    vector<pair<string, uint>> v;
    v.reserve(D.size());
    for (const auto & pair : D) {
        v.emplace_back(pair.first, pair.second);
    }
    sort(v.begin(), v.end(), cmp);

    vector<string> dict;
    vector<uint> remap = vector<uint>(P.size(), 0);
    for (uint i = 0; i < D.size(); i++) {
        dict.push_back(v[i].first);
        remap[v[i].second] = i;
    }

    vector<uint> parse = vector<uint>(P.size(), 0);
    for (uint i = 0; i < P.size(); i++) {
        parse[i] = remap[P[i]];
    }

    for (auto &it: dict) cout << it << endl;
    cout << parse << endl;

//    cout << parse << endl;
//    vector<uint> remap = {0, 5, 3, 6, 1, 2, 4, 0};
//    for (uint i = 0; i < parse.size(); i++) {
//        parse[i] = remap[parse[i]];
//    }
//    cout << parse << endl;
//
//    tfm_index<> tfm;
//    my_construct(tfm, parse);
//    cout << dot_repr_tfm(tfm) << endl;
//    wheeler_graph wg = wheeler_graph(tfm);
//    cout << wg.dot_repr() << endl;
//    wg_unparse(wg, dict);
//    cout << wg.dot_repr() << endl;
//    cout << reverse(wg_string(wg)) << endl;
//    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

