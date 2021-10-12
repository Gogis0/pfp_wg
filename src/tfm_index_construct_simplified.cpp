#include "functions.cpp"

class string_generator{
    vector<char> alphabet = {'A', 'C', 'G', 'T'};
    string s;

public:
    string_generator() { s += alphabet[alphabet.size() - 1]; }

    string next() {
        increase();
        return s;
    }

    string random(uint size) {
        randomize(size);
        return s;
    }

    void increase() {
        int inc = 1;
        for (int i = s.length() - 1; i >= 0; i--) {
            if (inc > 0) {
                uint pos = find(alphabet.begin(), alphabet.end(), s[i]) - alphabet.begin();
                if (pos == alphabet.size() - 1) {
                    s[i] = alphabet[0];
                    continue;
                } else {
                    pos++;
                    inc--;
                    s[i] = alphabet[pos];
                }
            }
        }
        if (inc > 0) {
            s = alphabet[0] + s;
        }
    }

    void randomize(uint size) {
        s = "";
        for (int i = 0; i < size; ++i)
            s += alphabet[rand() % alphabet.size()];
    }


};

int main() {
//    string_generator sg = string_generator();
//    for (uint i=0; i < 100; i++) {
//        string text = "$$" + sg.random(4500) + "##";
//        cout << text << endl;
//        vector<string> triggers = {"$$", "AC", "##"};
//
//        vector<string> dict;
//        vector<uint> parse;
//        fill_dict_and_parse(text, triggers, dict, parse);
//
//        tfm_index<> tfm = tfm_create(parse);
//        wheeler_graph wg = wheeler_graph(tfm);
//        wg_unparse(wg, dict);
//        position_end(wg, dict[0]);
//        wg_find_ordering(wg);
//
//        cout << (wg.is_valid()?"valid":"invalid") << endl;
//    }


    time_t timer = time(nullptr);
    stringstream message;

    // string path = "../BWT-Tunneling/testdata/pizzachili/random_bug1";
    // string T = "\x04\x04" + loadfile(path) + "\x03\x03";
//    string T = "\x04\x04" + string("AACCACACAAC") + "\x03\x03";
 //     vector<string> E = {"\x04\x04", "AC", "AG", "CG", "CT", "GT", "GA", "TA", "TC", "\x03\x03"};
//    vector<string> E = {"\x04\x04", "AC", "\x03\x03"};
 //    string T = "$$abeaaacda##";
 //    vector<string> E = {"$$", "aa", "##"};
    string T = "abeacdabeacda";
    vector<string> E = {"a"};

    vector<string> dict;
    vector<uint> parse;
    fill_dict_and_parse(T, E, dict, parse);
    cout << difftime(time(nullptr), timer) <<  " File parsed. " << dict.size() << " unique words. " << parse.size() << " total words." << endl;
    timer = time(nullptr);
    cout << parse << endl;
    for (uint i = 0; i < dict.size(); i++)
        cout << dict[i] << "\n";
    cout << parse.size() << endl;

    tfm_index<> tfm = tfm_create(parse);
    cout << difftime(time(nullptr), timer) << " tfm created" << endl;
    timer = time(nullptr);
    cout << tfm_repr(tfm) << endl;

    wheeler_graph wg = wheeler_graph(tfm);
    cout << difftime(time(nullptr), timer) << " tfm to wg" << endl;
    timer = time(nullptr);
    cout << wg.dot_repr() << endl;

    wg_unparse(wg, dict);
    cout << wg.dot_repr() << endl;
    position_end(wg, dict[0]);
    wg_find_ordering(wg);
    cout << difftime(time(nullptr), timer) << " unparsed and ordered" << endl;
    timer = time(nullptr);
    cout << wg.dot_repr() << endl;
    cout << (wg.is_valid()?"valid":"invalid") << endl;

    tfm_index<> tfm2 = wg_to_tfm(wg);
    cout << difftime(time(nullptr), timer) << " wg to tfm" << endl;

    return 0;
}

