rm -f ../bin/FreeBSD/x86_64/igrow ../obj/FreeBSD/x86_64/*.o
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86_64/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86_64/operation.o ../src/operation.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86_64/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86_64/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m64 -static -O3 -DNDEBUG -o ../bin/FreeBSD/x86_64/igrow ../obj/FreeBSD/x86_64/thread_pool.o ../obj/FreeBSD/x86_64/operation.o ../obj/FreeBSD/x86_64/ligand.o ../obj/FreeBSD/x86_64/main.o -L/home/hjli/boost_1_49_0/lib/x86_64 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
