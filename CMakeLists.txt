cmake_minimum_required(VERSION 3.10)
project(tunneling)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_FLAGS "-Wall -Wextra ")

include_directories(include)
link_directories(lib)
# include_directories(BWT-Tunneling/seqana/include)

# libraries
add_library(malloc_count bigrepair/ctph/malloc_count.h bigrepair/ctph/malloc_count.c)
add_library(utils bigrepair/ctph/utils.h bigrepair/ctph/utils.c)
add_library(basics bigrepair/largeb_repair/basics.h bigrepair/largeb_repair/basics.c)
add_library(records bigrepair/largeb_repair/records.h bigrepair/largeb_repair/records.c)
add_library(array bigrepair/largeb_repair/array.h bigrepair/largeb_repair/array.c)
add_library(hash bigrepair/largeb_repair/hash.h bigrepair/largeb_repair/hash.c)
add_library(heap bigrepair/largeb_repair/heap.h bigrepair/largeb_repair/heap.c)
target_link_libraries(heap array)

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

# gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000
add_executable(procdic bigrepair/procdic.c)
target_compile_definitions(procdic PUBLIC Unique=0x78000000)
add_executable(postproc bigrepair/postproc.c)
target_compile_definitions(postproc PUBLIC Unique=0x78000000)


# $(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS
add_executable(newscanNT bigrepair/ctph/newscan.cpp)
target_link_libraries(newscanNT dl malloc_count utils)
target_compile_definitions(newscanNT PUBLIC NOTHREADS)

# gcc $(CFLAGS) -o irepair irepair.o array.o hash.o heap.o records.o basics.o
add_executable(irepair bigrepair/largeb_repair/irepair.c)
target_link_libraries(irepair array hash heap records basics)

add_executable(despair bigrepair/repair/despair.c)

add_executable(tfm_explore src/tfm_explore.cpp)
target_link_libraries(tfm_explore sdsl divsufsort divsufsort64)

add_executable(tfm_index_construct_simplified src/tfm_index_construct_simplified.cpp)
target_link_libraries(tfm_index_construct_simplified sdsl divsufsort divsufsort64)