cmake_minimum_required(VERSION 3.10)
project(tunneling)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wextra ")

include_directories(include)
link_directories(lib)
include_directories(BWT-Tunneling/seqana/include)

add_executable(tfm_index_construct BWT-Tunneling/seqana/lib/tfm_index_construct.cpp)
target_link_libraries(
        tfm_index_construct
        sdsl divsufsort divsufsort64
)

add_executable(tfm_index_invert_to_binary BWT-Tunneling/seqana/lib/tfm_index_invert_to_binary.cpp)
target_link_libraries(
        tfm_index_invert_to_binary
        sdsl divsufsort divsufsort64
)