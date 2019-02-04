import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

n_points = 100
dx = 50.0/float(n_points)

data = np.loadtxt("density.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0], data[4*n_points:5*n_points-1,1]



def f_gauss(x,t):
    rho_0 = 10.0
    rho_1 = 50.0
    x_0 = 10.0
    sig = 2.0
    v = 1.0
    x_0t = x_0 + v*t
    return (v*(rho_1 - rho_0)*(((-2.0)*(x-x_0t)/(sig**2))*np.exp(-((x-x_0t)**2)/(sig**2))))

def f_r_shock(t):
    GAMMA = 5.0/3.0
    E = 100000/((GAMMA-1.0)*1.0)
    rho_1 = 1.0
    return 1.17*(E*t*t/rho_1)**0.2

plt.subplot(311)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

#times = np.array([0.00,2.00278,4.00166,6.00055,8.00333])
times = np.array([0.00,0.0200597,0.0400175,0.0600421,0.0800528])
lines = f_r_shock(times)

for l in lines:
    x_0 = 25.0
    plt.axvline(x=x_0+l,color='black')
    plt.axvline(x=x_0-l,color='black')

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

plt.subplot(312)

plt.plot(x0,rho0)
plt.plot(x1,rho1)
plt.plot(x2,rho2)
plt.plot(x3,rho3)
plt.plot(x4,rho4)

plt.ylabel("pressure [N/m^2]")
#plt.xlim([20,30])
#plt.ylim([0.001,1.5])

############################## velocity plot ##############################

data = np.loadtxt("velocity.txt")
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0]+0.5*dx, data[4*n_points:5*n_points-1,1]
x5,rho5 = data[5*n_points:6*n_points-1,0]+0.5*dx, data[5*n_points:6*n_points-1,1]

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