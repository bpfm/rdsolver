import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

n_points = 20
dx = 50.0/float(n_points)

data = np.loadtxt("density.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0], data[4*n_points:5*n_points-1,1]

rho_0 = 10.0
rho_1 = 50.0
x_0 = 10.0
sig = 2.0
v = 1.0

def f(x,t):
    x_0t = x_0 + v*t
    return (v*(rho_1 - rho_0)*(((-2.0)*(x-x_0t)/(sig**2))*np.exp(-((x-x_0t)**2)/(sig**2))))

plt.subplot(311)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

#times = np.array([0.00,2.00278,4.00166,6.00055,8.00333])
times = np.array([0.00,2.00278,4.00166,6.00055,8.00333])
x_0t = x_0 + times
lines = x_0t

#for l in lines:
        #plt.axvline(x=l,color='black')

#plt.xlabel("x [m]")
plt.ylabel("density [kg/m^3]")
#plt.xlim([20,30])
#plt.ylim([-1.0,1001])


############################## pressure plot ##############################

data = np.loadtxt("pressure.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0], data[4*n_points:5*n_points-1,1]
x5,rho5 = data[5*n_points:6*n_points-1,0], data[5*n_points:6*n_points-1,1]

plt.subplot(312)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

#plt.xlabel("x [m]")
plt.ylabel("pressure [N/m^2]")
#plt.xlim([20,30])
#plt.ylim([-1.0,1001])

############################## velocity plot ##############################

data = np.loadtxt("velocity.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0], data[4*n_points:5*n_points-1,1]
x5,rho5 = data[5*n_points:6*n_points-1,0], data[5*n_points:6*n_points-1,1]

plt.subplot(313)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

plt.xlabel("x [m]")
plt.ylabel("velocity [m/s]")
#plt.xlim([20,30])
#plt.ylim([-20,40])

plt.show()