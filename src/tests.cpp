#include <criterion/criterion.h>
#include "functions.cpp"

Test(external, empty_test) {}

Test(core, tfm_construction_small_example) {
    tfm_index<> tfm;
    auto parse = init_parse({1, 2, 1, 2});
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
    tfm_index<> tfm;
    auto parse = init_parse({2, 1, 4, 5, 3, 2, 1, 4, 5});   // this is easypeasy => a -> 1, e -> 2, ...
    my_construct(tfm, parse);

    wt_blcd<> wt = construct_from_vector({5, 2, 3, 0, 5, 1, 4});
    for (uint i=0; i < tfm.L.size(); i++) { cr_expect(tfm.L[i] == wt[i]); }
    bit_vector dout = int_vector<1>({1, 1, 1, 0, 1, 1, 1, 1});
    cr_expect(dout == tfm.dout);
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 1, 1, 0, 1});
    cr_expect(din == tfm.din);
}

Test(core, we_can_construct_the_same_tfm_from_parse_and_custom_vectors) {
    tfm_index<> tfm;
    auto parse = init_parse({1, 2, 1, 2});
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
    auto parse = init_parse({1, 2, 1, 2});

    tfm_index<> tfm;
    my_construct(tfm, parse);
    wheeler_graph wg = wheeler_graph(tfm);

    cr_assert(wg.dot_repr() == dot_repr_tfm(tfm));
}
