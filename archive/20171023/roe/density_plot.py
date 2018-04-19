import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

n_points = 200
dx = 50.0/float(n_points)

data = np.loadtxt("density.txt")
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0]+0.5*dx, data[4*n_points:5*n_points-1,1]

#print density

#x = np.zeros(shape = 50)
#y = np.zeros(shape = 50)

#for i in range(0,49):
#  x[i] = density[0][i]
#  y[i] = density[1][i]

plt.subplot(311)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

plt.xlabel("x [m]")
plt.ylabel("density [kg/m^3]")
#plt.xlim([20,30])
#plt.ylim([0.001,1.5])


############################## pressure plot ##############################


data = np.loadtxt("pressure.txt")
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0]+0.5*dx, data[4*n_points:5*n_points-1,1]
x5,rho5 = data[5*n_points:6*n_points-1,0]+0.5*dx, data[5*n_points:6*n_points-1,1]

#print density

#x = np.zeros(shape = 50)
#y = np.zeros(shape = 50)

#for i in range(0,49):
#  x[i] = density[0][i]
#  y[i] = density[1][i]

plt.subplot(312)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

plt.xlabel("x [m]")
plt.ylabel("pressure [N/m^2]")
#plt.xlim([20,30])
#plt.ylim([0.001,1.5])

#plt.show()

############################## velocity plot ##############################


data = np.loadtxt("velocity.txt")
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0]+0.5*dx, data[4*n_points:5*n_points-1,1]
x5,rho5 = data[5*n_points:6*n_points-1,0]+0.5*dx, data[5*n_points:6*n_points-1,1]

#print density

#x = np.zeros(shape = 50)
#y = np.zeros(shape = 50)

#for i in range(0,49):
#  x[i] = density[0][i]
#  y[i] = density[1][i]

plt.subplot(313)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

plt.xlabel("x [m]")
plt.ylabel("velocity [m/s]")
#plt.xlim([20,30])
#plt.ylim([0.001,1.5])

plt.show()