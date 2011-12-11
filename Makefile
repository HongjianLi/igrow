CC = g++ -O3 -DNDEBUG -std=gnu++0x

ifeq ($(TOOLSET), clang)
  CC = clang++ -static -O3 -DNDEBUG -std=gnu++11
else ifeq ($(TOOLSET), intel-linux)
  CC = icpc -static -O3 -DNDEBUG -std=gnu++0x
endif

bin/igrow: obj/atom.o obj/interact.o obj/bondlibrary.o obj/ligand.o obj/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

obj/%.o : src/%.cpp 
	$(CC) -o $@ -c $<

clean:
	rm -f bin/igrow obj/*.o
