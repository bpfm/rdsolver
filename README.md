# Residual Distribution Hydrodynamics Solver for Galaxy Formation Simulations

An implementation of the residual distriution (RD) partial differential equation (PDE) solver for the inviscid Eulerian fluid equations. 1D Roe hydro solver included for comparison with simple problems. We also include a standalone module for contructing periodic Delaunay meshes in both 2D and 3D, using the widely used CGAL library. A full description of the code can be found here (MNRAS paper): https://arxiv.org/pdf/2204.01757

In summary, `rdsolver` will evolve a compressible baryonic gas forward in time, within a static, periodic box, from some set of initial conditions (ICs). The intitial state requires the  density, velocity and internal energy distributions be defined across the box. At the time ICs are set up within the code, depending on the choice of case (outlined below). The domain is discretised into a static periodic Delaunay triangulation. The dual tesselation of this triangulation can be used to produce the piecewise uniform reconstruction of the fluid state at all positions, for the purposes of analyisis and visualisation.

Quickstart Guide

Compile Time Options

