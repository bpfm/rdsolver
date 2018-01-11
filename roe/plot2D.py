import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

data = np.loadtxt("density_slice.txt")

N_SIDE = 100
N_POINTS = N_SIDE*N_SIDE

x = data[0:N_POINTS,0]
y = data[0:N_POINTS,1]
rho = data[0:N_POINTS,2]

#triang = mtri.Triangulation(x,y)

#print triang

#plt.triplot(triang)
plt.scatter(x,y,rho)
plt.show()