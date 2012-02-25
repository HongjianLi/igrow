rm -f ../bin/FreeBSD/x86/igrow ../obj/FreeBSD/x86/*.o
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86/operation.o ../src/operation.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/FreeBSD/x86/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
clang++ -m32 -static -O3 -DNDEBUG -o ../bin/FreeBSD/x86/igrow ../obj/FreeBSD/x86/thread_pool.o ../obj/FreeBSD/x86/operation.o ../obj/FreeBSD/x86/ligand.o ../obj/FreeBSD/x86/main.o -L/home/hjli/boost_1_49_0/lib/x86 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
