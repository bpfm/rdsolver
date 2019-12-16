#g++ main.cpp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a
#./a.out

clang++ -Xpreprocessor -fopenmp -L /usr/local/opt/llvm/lib -lomp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a main.cpp -o lairds
./lairds


# cd output
# pythonw plot1D.py &

# cd output
# pythonw vor_plot2D.py

# cd ../../../exact_sod
# python sod_plot_rd.py &
