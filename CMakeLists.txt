cmake_minimum_required(VERSION 3.16)
project(tunneling)

set(CMAKE_CXX_STANDARD 11)

include_directories(~/projects/cds2/include)
link_directories(~/projects/cds2/lib)
include_directories(~/projects/cds_project/BWT-Tunneling/seqana/include)

add_executable(tfm_index_construct BWT-Tunneling/seqana/lib/tfm_index_construct.cpp)

target_link_libraries(
        tfm_index_construct
        sdsl divsufsort divsufsort64
)
