import matplotlib.pyplot as plt
import numpy as np

# verticle lines

xl = np.array([7.24,7.24])
yl = [0.0,4.0]

xl = xl + 25.0

n_points = 200
dx = 50.0/float(n_points)

data = np.loadtxt("density.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[9*n_points:10*n_points-1,0], data[9*n_points:10*n_points-1,1]

rho_0 = 10.0
rho_1 = 50.0
x_0 = 10.0
sig = 2.0
v = 1.0

def f(x,t):
    x_0t = x_0 + v*t
    return (v*(rho_1 - rho_0)*(((-2.0)*(x-x_0t)/(sig**2))*np.exp(-((x-x_0t)**2)/(sig**2))))

plt.subplot(311)

x_corr = x0#- 15.0

plt.plot(xl,yl,color = "black")

plt.plot(x_corr,rho0)
#plt.plot(x1,rho1)
#plt.plot(x2,rho2)
plt.plot(x_corr,rho3)
plt.plot(x_corr,rho4)

#times = np.array([0.00,2.00278,4.00166,6.00055,8.00333])
times = np.array([0.00,2.00278,4.00166,6.00055,8.00333])
x_0t = x_0 + times
lines = x_0t

#for l in lines:
        #plt.axvline(x=l,color='black')

#plt.xlabel("x [m]")
plt.ylabel("Density [kg/m^3]")
#plt.xlim([0,30])
plt.ylim([50-0.0051,50.0051])


############################## pressure plot ##############################

data = np.loadtxt("pressure.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
#x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
#x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[9*n_points:10*n_points-1,0], data[9*n_points:10*n_points-1,1]

plt.subplot(312)

plt.plot(x_corr,rho0)
#plt.plot(x1,rho1)
#plt.plot(x2,rho2)
plt.plot(x_corr,rho3)
plt.plot(x_corr,rho4)

#plt.xlabel("x [m]")
plt.ylabel("Pressure [N/m^2]")
#plt.xlim([0,30])
#plt.ylim([-1.0,1001])

############################## velocity plot ##############################

data = np.loadtxt("velocity.txt")
x0,rho0 = data[0:n_points-1,0], data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0], data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0], data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0], data[3*n_points:4*n_points-1,1]
x4,rho4 = data[9*n_points:10*n_points-1,0], data[9*n_points:10*n_points-1,1]

plt.subplot(313)

plt.plot(x_corr,rho0)
#plt.plot(x1,rho1)
#plt.plot(x2,rho2)
plt.plot(x_corr,rho3)
plt.plot(x_corr,rho4)

plt.xlabel("x [m]")
plt.ylabel("X Velocity [m/s]")
#plt.xlim([0,30])
#plt.ylim([0,20])

plt.show()
