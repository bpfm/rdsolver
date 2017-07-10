import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

data = np.loadtxt("density.txt")                       # myfile.txt contains 4 columns of numbers
x0,rho0 = data[0:99,0], data[0:99,1]
x1,rho1 = data[100:199,0], data[100:199,1]
x2,rho2 = data[200:299,0], data[200:299,1]
x3,rho3 = data[300:399,0], data[300:399,1]

#print density

#x = np.zeros(shape = 50)
#y = np.zeros(shape = 50)

#for i in range(0,49):
#	x[i] = density[0][i]
#	y[i] = density[1][i]

plt.semilogy(x0,rho0)
plt.semilogy(x1,rho1)
plt.semilogy(x2,rho2)
plt.semilogy(x3,rho3)
plt.show()