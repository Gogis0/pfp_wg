#include "functions.cpp"

int main() {
    tfm_index<> tfm;
    // auto parse = init_parse({1, 2, 1, 2});
    auto parse = init_parse({2, 1, 4, 5, 3, 2, 1, 4, 5});
    my_construct(tfm, parse);
    print_tfm(tfm);
    print_original(tfm);

    return 0;
}

