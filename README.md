# Residual Distribution Hydrodynamics Solver for Galaxy Formation Simulations

An implementation of the residual distriution (RD) partial differential equation (PDE) solver for the inviscid Eulerian fluid equations. Included as an isolated tool is code to generate periodic Delaunay triangulations using the CGAL mesh library.

Residual distribution solver for 2D Euler equations:
- 1st order LDA 
- 1st and 2nd order N
- 1st and 2nd order B

Delaunay triangulation construction routines (seperate from main code for now) using
- CGAL (periodic if required)


