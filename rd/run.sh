g++ main.cpp -I OpenBLAS/include/ OpenBLAS/lib/libopenblas.a
./a.out

cd output
python plot1D.py &

# cd ../../../exact_sod
# python sod_plot_rd.py &
