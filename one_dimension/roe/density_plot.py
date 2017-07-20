import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

n_points = 50
dx = 20.0/float(n_points)

data = np.loadtxt("density.txt")                       # myfile.txt contains 4 columns of numbers
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]

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
plt.xlabel("x [m]")
plt.ylabel("density [kg/m^3]")
plt.xlim([8,12])
plt.show()