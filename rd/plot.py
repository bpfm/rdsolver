import numpy as np
import matplotlib.pylab as plt
import matplotlib.tri as mtri

data = np.loadtxt("positions.txt")

x = data[:,0]
y = data[:,1]

triang = mtri.Triangulation(x,y)

print triang

plt.triplot(triang)
plt.show()