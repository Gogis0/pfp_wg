#include "functions.cpp"

int main() {
    tfm_index<> tfm;
    // auto parse = init_parse({1, 2, 1, 2});
    auto parse = init_parse({1, 2, 1, 3, 1, 2, 1, 3, 1});
    my_construct(tfm, parse);
    print_tfm(tfm);
    print_original(tfm);

    tfm_index<> tfm2;
    wt_blcd<> wt = construct_from_vector({2, 2, 0, 1, 1, 0, 0, 0});
    bit_vector dout = int_vector<1>({1, 1, 1, 1, 0, 1, 1, 1, 1});
    bit_vector din  = int_vector<1>({1, 1, 1, 1, 1, 1, 1, 0, 1});
    uint64_t text_len = 9;
    construct_tfm_index_tmp(tfm2, text_len, move(wt), move(dout), move(din));
    print_tfm(tfm2);
    print_original(tfm2);

    return 0;
}

