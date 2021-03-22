#g++ main.cpp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a
#./a.out

# clang++ -Xpreprocessor -fopenmp -L /usr/local/opt/llvm/lib -lomp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a main.cpp -o lairds

clang++ -Xpreprocessor -fopenmp -lomp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a main.cpp -o lairds
./lairds

# clang++ -I /home/morton/local/include  -L /home/morton/local/lib -lopenblas -o lairds main.cpp

