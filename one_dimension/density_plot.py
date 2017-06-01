import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

data = np.loadtxt("density.txt")                       # myfile.txt contains 4 columns of numbers
x0,rho0 = data[0:49,0], data[0:49,1]
x1,rho1 = data[51:100,0], data[51:100,1]
x2,rho2 = data[102:151,0], data[102:151,1]

#print density

#x = np.zeros(shape = 50)
#y = np.zeros(shape = 50)

#for i in range(0,49):
#	x[i] = density[0][i]
#	y[i] = density[1][i]

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.show()