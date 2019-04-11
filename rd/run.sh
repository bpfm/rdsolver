g++ main.cpp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a
./a.out

cd output
pythonw plot1D.py &

# cd output
# pythonw vor_plot2D.py

# cd ../../../exact_sod
# python sod_plot_rd.py &
