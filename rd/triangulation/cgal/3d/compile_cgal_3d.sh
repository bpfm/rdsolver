PATH_TO_CGAL_SCRIPTS="/home/morton/CGAL-5.0.2/Scripts/scripts/cgal_create_CMakeLists"
bash $PATH_TO_CGAL_SCRIPTS -s cgal_periodic3D
cmake -DCMAKE_BUILD_TYPE=Release .
