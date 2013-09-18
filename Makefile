CC=clang++ -std=c++11 -O3

bin/igrow: obj/thread_pool.o obj/operation.o obj/ligand.o obj/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

obj/%.o: src/%.cpp
	$(CC) -o $@ $< -c

clean:
	rm -f bin/igrow obj/*.o
