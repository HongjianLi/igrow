BOOST_ROOT = $(HOME)/boost_1_52_0
CC = g++ -O3 -DNDEBUG -std=c++11

ifeq ($(TOOLSET), clang)
  CC = clang++ -O3 -DNDEBUG -std=c++11
else ifeq ($(TOOLSET), intel)
  CC = icpc -O3 -DNDEBUG -std=c++11
endif

bin/igrow: obj/thread_pool.o obj/operation.o obj/ligand.o obj/main.o
	$(CC) -o $@ $^ -L$(BOOST_ROOT)/lib/x86_64 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -I$(BOOST_ROOT) -c

clean:
	rm -f bin/igrow obj/*.o
