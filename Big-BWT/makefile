# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

# main executables 
EXECS=bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pscan.x pfbwt.x pfbwt64.x unparse remap
# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscanNT.x pfbwtNT.x pfbwtNT64.x

# targets not producing a file declared phony
.PHONY: all clean tarfile

all: $(EXECS) $(EXECS_NT)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

bwtparse: bwtparse.c gsa/gsacak.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

bwtparse64: bwtparse.c gsa/gsacak64.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl -DM64

simplebwt: simplebwt.c gsa/gsacak.o
	$(CC) $(CFLAGS) -o $@ $^

simplebwt64: simplebwt.c gsa/gsacak64.o
	$(CC) $(CFLAGS) -o $@ $^ -DM64

# cnewscan executable to scan gzipped files (currently not active) 
cnewscan.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl -DNOTHREADS

newscanNT.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

newscan.x: newscan.cpp newscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan.cpp malloc_count.o utils.o xerrors.o -ldl -lz -pthread

pscan.x: pscan.cpp pscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ pscan.cpp malloc_count.o utils.o xerrors.o -ldl -pthread


# prefix free BWT construction
pfbwt.x: pfbwt.cpp pfthreads.hpp gsa/gsacak.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak.o utils.o xerrors.o malloc_count.o -pthread -ldl

pfbwt64.x: pfbwt.cpp pfthreads.hpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o -pthread -ldl -DM64

# TO BE REMOVED? (now pfbwt*.x works with malloc_count)
# prefix free BWT construction without threads: useful since supports malloc_count
pfbwtNT.x: pfbwt.cpp gsa/gsacak.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

# as above for large files
pfbwtNT64.x: pfbwt.cpp gsa/gsacak64.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS -DM64

unparse: unparse.c utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

remap: remap.c
	$(CC) $(CFLAGS) -o $@ $< -O -g -Wall -lm 

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.[ch]pp pscan.[ch]pp pfbwt.cpp pfthreads.hpp simplebwt.c bwtparse.c unparse.c remap.c makefile utils.[ch] xerrors.[ch] f2s.py gsa/gsacak.[ch] gsa/LICENSE gsa/README.md malloc_count.[ch]

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
	rm -f *.sdsl
	rm -f *.x

include ~/sdsl-lite/Make.helper
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

#DBGEXECUTABLES = tfm_index_construct.x tfm_index_invert.x tfm_index_invert_to_binary.x dbg_edgespectrum.x
DBGEXECUTABLES =  tfm_index_construct.x tfm_index_invert.x tfm_index_invert_to_binary.x tfm_generate_aux.x
#TRIEEXECUTABLES = create_trie_input.x trie_construct.x mp_search.x
TRIEEXECUTABLES = create_trie_input.x  mp_search.x

all: $(DBGEXECUTABLES) $(TRIEEXECUTABLES)

#DBG EXECUTABLES
tfm_index_construct.x: include/tfm_index.hpp include/dbg_algorithms.hpp tfm_index_construct.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib tfm_index_construct.cpp -o tfm_index_construct.x $(LIBS)

tfm_index_invert.x: include/tfm_index.hpp include/dbg_algorithms.hpp tfm_index_invert.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib tfm_index_invert.cpp -o tfm_index_invert.x $(LIBS)

tfm_generate_aux.x: include/tfm_index.hpp tfm_generate_aux.cpp utils.o
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib utils.o tfm_generate_aux.cpp -o tfm_generate_aux.x $(LIBS)

tfm_index_invert_to_binary.x: include/tfm_index.hpp include/dbg_algorithms.hpp tfm_index_invert_to_binary.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib tfm_index_invert_to_binary.cpp -o tfm_index_invert_to_binary.x $(LIBS)

# dbg_edgespectrum.x does not work yet due to char->int alphabet transformation
#dbg_edgespectrum.x: include/dbg_algorithms.hpp lib/dbg_edgespectrum.cpp
#	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
#	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib lib/dbg_edgespectrum.cpp -o dbg_edgespectrum.x $(LIBS)

#TRIE EXECUTABLES
create_trie_input.x: create_trie_input.cpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib create_trie_input.cpp -o create_trie_input.x $(LIBS)

trie_construct.x: trie_construct.cpp include/tries.hpp include/trie_xbwt.hpp include/trie_txbwt.hpp include/trie_construct_algorithms.hpp include/dbg_algorithms.hpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib trie_construct.cpp -o trie_construct.x $(LIBS)

mp_search.x: mp_search.cpp include/tries.hpp include/trie_xbwt.hpp include/trie_txbwt.hpp
	$(MY_CXX) -Wall -Wextra $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(C_OPTIONS) \
	-I$(INC_DIR) -L$(LIB_DIR) -Iinclude -Llib mp_search.cpp -o mp_search.x $(LIBS)

