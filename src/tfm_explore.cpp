#include "../BWT-Tunneling/seqana/include/tfm_index.hpp"

using namespace std;

int main() {
    tfm_index<> tfm;
    load_from_file(tfm, "../data/yeast.fasta.tunnel");

    const tfm_index<>::wt_type &L = tfm.L;
    const std::vector<sdsl::int_vector<>::size_type> &C = tfm.C;
    const tfm_index<>::bit_vector_type &dout = tfm.dout;
//    const rank_type &dout_rank = m_dout_rank;
//    const select_type &dout_select = m_dout_select;
    const tfm_index<>::bit_vector_type &din = tfm.din;
//    const rank_type &din_rank = m_din_rank;
//    const select_type &din_select = m_din_select;

    tfm_index<>::size_type size = tfm.size();
    tfm_index<>::nav_type nav = tfm.end();
    tfm_index<>::value_type value = tfm.preceding_char(nav);
//    value_type backwardstep(nav_type &pos)
//    size_type serialize(
//            std::ostream &out, sdsl::structure_tree_node *v, std::string name
//    )
//    void load(std::istream &in)
    cout << L.size() << "\n";
    cout << L.sigma << "\n";
//    for (uint i = 0; i < L.bv.size(); i++)
//        cout << L[i] << " ";
//    cout << endl;
}

