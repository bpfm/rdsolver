import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

data = np.loadtxt("positions.txt")

x = data[:,0]
y = data[:,1]
rho = data[:,2]

#triang = mtri.Triangulation(x,y)

print triang

#plt.triplot(triang)
plt.scatter(x,y,rho)
plt.show()