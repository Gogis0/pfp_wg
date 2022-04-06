# compilation flags
CC=/usr/bin/cc
CFLAGS=-O3 -Wall -std=c99 -g

CXX=/usr/bin/c++
CXX_FLAGS=-std=c++11 -Wall -Wextra -DNDEBUG
CXX_OPT_FLAGS=-O3 -ffast-math -funroll-loops -msse4.2 -march=native -DHAVE_CXA_DEMANGLE

LIBS = -lsdsl -ldivsufsort -ldivsufsort64

# main executables
EXECS=pfwg.x newscanNT.x tfm_index_construct.x tfm_index_invert.x

TEST_INPUTS = data/yeast.raw data/yeast.fasta
INVERT_TESTS = test/inv.00 test/inv.01 test/inv.02 test/inv.03 test/inv.04 test/inv.05

# targets not producing a file declared phony
.PHONY: all clean tarfile autotest

all: $(EXECS) $(EXECS_NT) $(DBGEXECUTABLES) $(TRIEEXECUTABLES)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

newscanNT.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

# prefix free WG construction
pfwg.x: pfwg.cpp include/tfm_index.hpp  gsa/gsacak.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfwg.cpp -Iinclude -Llib gsa/gsacak.o  utils.o xerrors.o malloc_count.o -pthread -ldl $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

#DBG EXECUTABLES
tfm_index_construct.x: include/tfm_index.hpp include/dbg_algorithms.hpp tfm_index_construct.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_OPT_FLAGS) $(C_OPTIONS) \
	-Iinclude -Llib tfm_index_construct.cpp -o tfm_index_construct.x $(LIBS)

tfm_index_invert.x: include/tfm_index.hpp include/dbg_algorithms.hpp tfm_index_invert.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_OPT_FLAGS) $(C_OPTIONS) \
	-Iinclude -Llib tfm_index_invert.cpp -o tfm_index_invert.x $(LIBS)

autotest: $(INVERT_TESTS)
	@echo "\033[0;32mAll tests passed!\033[0m"

$(TEST_INPUTS): tfm_index_invert.x $(EXECS) 
	./bigbwt $@
	./tfm_index_invert.x $@.tunnel
	diff $@.parse $@.untunneled

$(INVERT_TESTS): tfm_index_invert.x
	./tfm_index_invert.x $@
	diff $@ $@.untunneled

num_runs: num_runs.c
	gcc num_runs.c -o num_runs

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
	rm -f *.x

