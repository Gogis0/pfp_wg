#include <criterion/criterion.h>
#include "functions.cpp"

Test(external, empty_test) {}

Test(core, tfm_construction_small_example) {
    vector<uint> parse = {1, 2, 1, 2};

    tfm_index<> tfm;
    my_construct(tfm, parse);

    wt_blcd<> wt = construct_from_vector({2, 2, 0, 1});
    for (uint i=0; i < tfm.L.size(); i++) {
        cr_expect(tfm.L[i] == wt[i]);
    }

    bit_vector dout = int_vector<1>({1, 1, 0, 1, 1});
    cr_expect(dout == tfm.dout);
    bit_vector din  = int_vector<1>({1, 1, 1, 0, 1});
    cr_expect(din == tfm.din);
}

Test(core, tfm_construction_easypeasy) {
    // this is            e  a  s  y  p  e  a  s  y
    vector<uint> parse = {2, 1, 4, 5, 3, 2, 1, 4, 5};
    tfm_index<> tfm;
    my_construct(tfm, parse);

    wt_blcd<> wt = construct_from_vector({5, 2, 3, 0, 5, 1, 4});
    for (uint i=0; i < tfm.L.size(); i++) { cr_expect(tfm.L[i] == wt[i]); }
    bit_vector dout = int_vector<1>({1, 1, 1, 0, 1, 1, 1, 1});
    cr_expect(dout == tfm.dout);
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 1, 1, 0, 1});
    cr_expect(din == tfm.din);
}

Test(core, we_can_construct_the_same_tfm_from_parse_and_custom_vectors) {
    vector<uint> parse = {1, 2, 1, 2};

    tfm_index<> tfm;
    my_construct(tfm, parse);

    tfm_index<> tfm2;
    wt_blcd<> wt = construct_from_vector({2, 2, 0, 1});
    bit_vector dout = int_vector<1>({1, 1, 0, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 0, 1});
    uint64_t text_len = tfm.size();

    construct_tfm_index_tmp(tfm2, text_len, move(wt), move(dout), move(din));

    for (uint i=0; i < tfm.L.size(); i++) { cr_expect(tfm.L[i] == tfm2.L[i]); }
    cr_expect(tfm2.dout == tfm.dout);
    cr_expect(tfm2.din == tfm.din);
}

Test(core, naive_representation_is_the_same_as_tfm_representation) {
    vector<uint> parse = {1, 2, 1, 2};

    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);

    cr_assert(wg.dot_repr() == dot_repr_tfm(tfm));
}

Test(core, we_can_read_dictionary) {
    vector<string> dict = read_dict("../data/yeast.fasta.dict");
    // for (const string &s: dict) {cout << s << endl;}
    cr_assert(dict[0] == "\u0002CCACACC");
    cr_assert(dict[1] == "CACCACACACC");
    cr_assert(dict[2] == "CACCACACC");
    cr_assert(dict[3] == "CACCCACACACACA\u0002\u0002\u0002\u0002");
    cr_assert(dict[4] == "CACCCACACACACACC");
    cr_assert(dict[5] == "CACCCACACACC");
}

Test(core, relabel_works) {
    wheeler_graph wg = wheeler_graph(2);
    wg.add_edge(0, 1, 0, 0);
    vector<string> dict = {"ABAABAB"};
    expand_edge(wg, 0, dict[0]);
    cr_assert(
        "digraph G {\n"
        "\t0 -> 2 [label = \"t=0\\nl=65\"]\n"
        "\t2 -> 3 [label = \"t=0\\nl=66\"]\n"
        "\t3 -> 4 [label = \"t=0\\nl=65\"]\n"
        "\t4 -> 5 [label = \"t=0\\nl=65\"]\n"
        "\t5 -> 6 [label = \"t=0\\nl=66\"]\n"
        "\t6 -> 7 [label = \"t=0\\nl=65\"]\n"
        "\t7 -> 1 [label = \"t=0\\nl=66\"]\n"
        "}" == wg.dot_repr()
    );
}

Test(core, expand_graph_returns_sensible_answer) {
    vector<uint> parse = {1, 2, 1, 2};  // aba aca aba aca
    vector<string> dict = {"a", "ab", "ac"};

    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    cr_assert(
        "digraph G {\n"
        "\t0 -> 4 [label = \"t=0\\nl=97\"]\n"
        "\t4 -> 3 [label = \"t=0\\nl=99\"]\n"
        "\t3 -> 5 [label = \"t=0\\nl=97\"]\n"
        "\t5 -> 1 [label = \"t=1\\nl=98\"]\n"
        "\t1 -> 6 [label = \"t=1\\nl=97\"]\n"
        "\t6 -> 3 [label = \"t=0\\nl=99\"]\n"
        "\t3 -> 7 [label = \"t=0\\nl=97\"]\n"
        "\t7 -> 2 [label = \"t=0\\nl=98\"]\n"
        "\t2 -> 0 [label = \"t=0\\nl=97\"]\n"
        "}" == wg.dot_repr()
    );
}

Test(core, finds_correct_ordering) {
    vector<uint> parse = {1, 2, 1, 2};  // aba aca aba aca
    vector<string> dict = {"a", "ab", "ac"};
    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    wg_unparse(wg, dict);
    wg_find_ordering(wg);
    cr_assert(wg.is_valid());
}

Test(core, test_our_end) {
    vector<uint> parse = {1, 2, 1, 2};
    tfm_index<> tfm;
    my_construct(tfm, parse);
    auto p = tfm.our_end();
    for (uint i = 0; i < tfm.size(); i++) {
        // cout << tfm.preceding_char(p) << endl;
        tfm.backwardstep(p);
    }
}

Test(core, test_backward) {
    vector<uint> parse = {1, 2, 1, 2};
    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    auto p = wg.end();
    cr_assert(p.first == 0 && p.second == 0);
    wg.backward(p);
    cr_assert(p.first == 3 && p.second == 0);
    wg.backward(p);
    cr_assert(p.first == 1 && p.second == 0);
    wg.backward(p);
    cr_assert(p.first == 3 && p.second == 1);
    wg.backward(p);
    cr_assert(p.first == 2 && p.second == 0);
    wg.backward(p);
    cr_assert(p.first == 0 && p.second == 0);
}

Test(core, test_forward) {
    vector<uint> parse = {1, 2, 1, 2};
    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);
    auto p = wg.end();
    cr_assert(p.first == 0 && p.second == 0);
    wg.forward(p);
    cr_assert(p.first == 2 && p.second == 0);
    wg.forward(p);
    cr_assert(p.first == 3 && p.second == 1);
    wg.forward(p);
    cr_assert(p.first == 1 && p.second == 0);
    wg.forward(p);
    cr_assert(p.first == 3 && p.second == 0);
    wg.forward(p);
    cr_assert(p.first == 0 && p.second == 0);
}
