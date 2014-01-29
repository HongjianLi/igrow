CC=clang++ -std=c++11 -O2

bin/igrow: obj/io_service_pool.o obj/operation.o obj/ligand.o obj/safe_counter.o obj/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_filesystem -lboost_program_options

obj/%.o: src/%.cpp
	$(CC) -o $@ $< -c

clean:
	rm -f bin/igrow obj/*.o
