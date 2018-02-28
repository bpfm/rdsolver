import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
N_SIDE = 50
N_POINTS = N_SIDE*N_SIDE

SLIDE = 0

# read in density values
points0 = np.zeros(shape=[N_POINTS, 2])
points1 = np.zeros(shape=[N_POINTS, 2])

data = np.loadtxt("density.txt")

points0[:,0] = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,0]
points0[:,1] = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,1]
rho0 = data[SLIDE*N_POINTS:(SLIDE+1)*N_POINTS,2]

points1[:,0] = data[(SLIDE+1)*N_POINTS:(SLIDE+2)*N_POINTS,0]
points1[:,1] = data[(SLIDE+1)*N_POINTS:(SLIDE+2)*N_POINTS,1]
rho1 = data[(SLIDE+1)*N_POINTS:(SLIDE+2)*N_POINTS,2]

change = rho1 - rho0

# generate Voronoi tessellation
vor = Voronoi(points0)

# normalize chosen colormap
norm = mpl.colors.Normalize(vmin=np.amin(change), vmax=np.amax(change), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Oranges)

# plot Voronoi diagram, and fill finite regions with color mapped from speed value
voronoi_plot_2d(vor,show_points=False,show_vertices=False,line_width=0)
for r in range(len(vor.point_region)):
    region = vor.regions[vor.point_region[r]]
    if not -1 in region:
        polygon = [vor.vertices[i] for i in region]
        plt.fill(*zip(*polygon), color=mapper.to_rgba(change[r]))

plt.xlabel("x [m]")
plt.ylabel("y [m]")
#plt.colorbar()

plt.show()
