cmake_minimum_required(VERSION 3.10)
project(tunneling)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)
link_directories(lib)
include_directories(BWT-Tunneling/seqana/include)

add_executable(tfm_index_construct BWT-Tunneling/seqana/lib/tfm_index_construct.cpp)

target_link_libraries(
        tfm_index_construct
        sdsl divsufsort divsufsort64
)
