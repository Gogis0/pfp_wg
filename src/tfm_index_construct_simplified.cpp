#include "functions.cpp"

int main() {
    time_t timer = time(nullptr);
    stringstream message;

    string path = "../BWT-Tunneling/testdata/pizzachili/dna.small";
    string T = "\x04\x04" + loadfile(path) + "\x03\x03";
    // vector<string> E = {"\x04\x04", "AC", "AG", "CG", "CT", "GT", "GA", "TA", "TC", "\x03\x03"};
    vector<string> E = {"\x04\x04", "AC", "\x03\x03"};

    vector<string> dict;
    vector<uint> parse;
    fill_dict_and_parse(T, E, dict, parse);
    cout << difftime(time(nullptr), timer) <<  " File parsed. " << dict.size() << " unique words. " << parse.size() << " total words." << endl;
    timer = time(nullptr);

    tfm_index<> tfm = tfm_create(parse);
    cout << difftime(time(nullptr), timer) << " tfm created" << endl;
    timer = time(nullptr);

    wheeler_graph wg = wheeler_graph(tfm);
    cout << difftime(time(nullptr), timer) << " tfm to wg" << endl;
    timer = time(nullptr);

    wg_unparse(wg, dict);
    position_end(wg, dict[0]);
    wg_find_ordering(wg);
    cout << difftime(time(nullptr), timer) << " unparsed and ordered" << endl;
    timer = time(nullptr);

    tfm_index<> tfm2 = wg_to_tfm(wg);
    cout << difftime(time(nullptr), timer) << " wg to tfm" << endl;

    return 0;
}

