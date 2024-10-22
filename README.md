# Residual Distribution Hydrodynamics Solver for Galaxy Formation Simulations

An implementation of the residual distriution (RD) partial differential equation (PDE) solver for the inviscid Eulerian fluid equations. 1D Roe hydro solver included for comparison with simple problems. We also include a standalone module for contructing periodic Delaunay meshes in both 2D and 3D, using the widely used CGAL library. A full description of the code can be found here (MNRAS paper): https://arxiv.org/pdf/2204.01757

In summary, `rdsolver` will evolve a compressible baryonic gas forward in time, within a static, periodic box, from some set of initial conditions (ICs). The intitial state requires the  density, velocity and internal energy distributions be defined across the box. At the time ICs are set up within the code, depending on the choice of case (outlined below). The domain is discretised into a static periodic Delaunay triangulation. The dual tesselation of this triangulation can be used to produce the piecewise uniform reconstruction of the fluid state at all positions, for the purposes of analyisis and visualisation.

## Prerequisits

OpenBLAS

OpenMP

C++ Compiler (tested most widely with clang)

CGAL (stand alone to generate meshes)


## Run Guide

1. Download repository using `git clone https://github.com/bpfm/rdsolver.git`
2. Move to rd directory
3. Edit CMakeLists.txt with the appropriate setting for your machine:
    -  Compiler
    -  OpenBLAS and OpenMP include paths
    -  OpenBLAS lib directory
4. Run `cmake .` to generate `Makefile`
5. Select 2D, 3D or both cases to prepare by commenting out first/second `add_executable` lines
6. Set compile time settings by copying constants.h.TEMPLATE or constants3D.h.TEMPLATE to create constants.h etc.
7. Edit to choose scenario (selecting case by uncommenting from first block, see below for details)
8. Generate mesh with appropriate dimensions
    - Move to triangulation/cgal/2d or triangulation/cgal/3d
    - Edit CMakeLists.txt with your CGAL location
    - Run `compile_cgal_2d.sh` or `compile_cgal_3d.sh` (whichever is appropriate to your chosen scenario)
    - Edit `cgal_periodic_2d.cpp` or `cgal_periodic_3d.cpp` to give it the correct edge length for the chosen problem, and the correct number of vertices
    - Run `run_cgal_2d.sh` or `run_cgal_3d.sh`
    - Copy output `.txt` file to `rd` directory, or whereever you plan to run
    - Move back to `rd`
9. Compile by running Makefile with `make`
10. Run executable (or move to run directory and run)
  
## Compile Time Options

