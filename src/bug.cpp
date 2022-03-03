#include "functions.cpp"

void print_parse(
    const string &text, const vector<string> &dict, const vector<uint> &parse
) { 
    cout << text << "\n";
    for (uint i = 0; i < parse.size(); i++) {
        printf("%-*d", (int)dict[parse[i]].length(), parse[i]);
    }
    cout << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <textfile>\n", argv[0]);
        exit(1);
    }

    string path = argv[1];
    string text = "$$" + loadfile(path) + "##";

    vector<string> triggers = {"$$", "AC", "##"};

    vector<string> dict;
    vector<uint> parse;
    fill_dict_and_parse(text, triggers, dict, parse);

    print_parse(text, dict, parse);

    tfm_index<> tfm = tfm_create(parse);
    cout << tfm_repr(tfm) << endl;

    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    position_end(wg, dict[0]);
    wg_find_ordering(wg);
    cout << wg.dot_repr() << endl;

    cout << (wg.is_valid()?"valid":"invalid") << endl;

    return 0;
}

