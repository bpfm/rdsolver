import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

N_SIDE = 200
N_POINTS = N_SIDE*N_SIDE

SLIDE = 4

# read in density values
points = np.zeros(shape=[N_POINTS, 2])

data = np.loadtxt("density.txt")

points[:,0] = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,0]
points[:,1] = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,1]
rho = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,2]

# generate Voronoi tessellation
vor = Voronoi(points)

# normalize chosen colormap
norm = mpl.colors.Normalize(vmin=np.amin(rho), vmax=np.amax(rho), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Oranges)

# plot Voronoi diagram, and fill finite regions with color mapped from speed value
voronoi_plot_2d(vor,show_points=False,show_vertices=False,line_width=0)
for r in range(len(vor.point_region)):
    region = vor.regions[vor.point_region[r]]
    if not -1 in region:
        polygon = [vor.vertices[i] for i in region]
        plt.fill(*zip(*polygon), color=mapper.to_rgba(rho[r]))

plt.xlabel("x [m]")
plt.ylabel("y [m]")

plt.show()
