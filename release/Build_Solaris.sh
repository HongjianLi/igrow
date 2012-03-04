rm -f ../bin/Solaris/x86_64/igrow ../obj/Solaris/x86_64/*.o
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/operation.o ../src/operation.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -o ../bin/Solaris/x86_64/igrow ../obj/Solaris/x86_64/thread_pool.o ../obj/Solaris/x86_64/operation.o ../obj/Solaris/x86_64/ligand.o ../obj/Solaris/x86_64/main.o -L/home/hjli/boost_1_49_0/lib/x86_64 -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
rm -f ../bin/Solaris/x86/igrow ../obj/Solaris/x86/*.o
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/operation.o ../src/operation.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -o ../bin/Solaris/x86/igrow ../obj/Solaris/x86/thread_pool.o ../obj/Solaris/x86/operation.o ../obj/Solaris/x86/ligand.o ../obj/Solaris/x86/main.o -L/home/hjli/boost_1_49_0/lib/x86 -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
