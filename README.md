# Residual Distribution Hydrodynamics Solver for Galaxy Formation Simulations

An implementation of the residual distriution (RD) partial differential equation (PDE) solver for the inviscid Eulerian fluid equations. 1D Roe hydro solver included for comparison with simple problems. Included as an isolated tool is code to generate periodic Delaunay triangulations using the CGAL mesh library. Methods and tests presented in MNRAS paper here: https://arxiv.org/pdf/2204.01757

In summary, `rdsolver` will evolve a compressible baryonic gas forward in time, within a static, periodic box, from some set of initial conditions. The intitial state requires the  density, velocity and internal energy distributions be defined across the box. The domain is discretised into a static Delaunay triangulation. The dual tesselation of this triangulation is used to produce the piecewise uniform reconstruction of the fluid state at all positions.

Residual distribution solver for 2D Euler equations:
- 1st order LDA 
- 1st and 2nd order N
- 1st and 2nd order B



