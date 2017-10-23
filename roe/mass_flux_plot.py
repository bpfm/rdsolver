import matplotlib.pyplot as plt
import numpy as np

rho_0 = 10.0
rho_1 = 50.0
x_0 = 10.0
sig = 2.0
v = 1.0
dt = 0.00974115 # 0.00389646 # 0.0194823 # 0.00974115

n_points = 200
dx = 50.0/float(n_points)

def f(x,t):
    x_0t = x_0 + v*t
    return (v*(rho_1 - rho_0)*(((-2.0)*(x-x_0t)/(sig**2))*np.exp(-((x-x_0t)**2)/(sig**2))))

data = np.loadtxt("density.txt")
x0,rho0 = data[0:n_points-1,0]+0.5*dx, data[0:n_points-1,1]
x1,rho1 = data[n_points:2*n_points-1,0]+0.5*dx, data[n_points:2*n_points-1,1]
x2,rho2 = data[2*n_points:3*n_points-1,0]+0.5*dx, data[2*n_points:3*n_points-1,1]
x3,rho3 = data[3*n_points:4*n_points-1,0]+0.5*dx, data[3*n_points:4*n_points-1,1]
x4,rho4 = data[4*n_points:5*n_points-1,0]+0.5*dx, data[4*n_points:5*n_points-1,1]

x5,rho5 = data[5*n_points:6*n_points-1,0]+0.5*dx, data[5*n_points:6*n_points-1,1]
x6,rho6 = data[6*n_points:7*n_points-1,0]+0.5*dx, data[6*n_points:7*n_points-1,1]
x7,rho7 = data[7*n_points:8*n_points-1,0]+0.5*dx, data[7*n_points:8*n_points-1,1]
x8,rho8 = data[8*n_points:9*n_points-1,0]+0.5*dx, data[8*n_points:9*n_points-1,1]
x9,rho9 = data[9*n_points:10*n_points-1,0]+0.5*dx, data[9*n_points:10*n_points-1,1]

xface0 = x0-dx

change = (rho0-rho1)/dt
change2 = (rho1-rho2)/dt
change3 = (rho2-rho3)/dt
change4 = (rho3-rho4)/dt

change5 = (rho4-rho5)/dt
change6 = (rho5-rho6)/dt
change7 = (rho6-rho7)/dt
change8 = (rho7-rho8)/dt
change9 = (rho8-rho9)/dt


plt.figure(1)
#plt.subplot(311)

#plt.xlim(0,20)
#plt.ylabel("density flow [kg/m^3/s]")

#plt.plot(x0,change)
#plt.plot(x0,f(xface0,0.0))


#plt.subplot(312)

#plt.xlim(0,20)
#plt.plot(x0,change)
#plt.plot(x0,f(xface0,dt))


#plt.subplot(313)

diff = (f(xface0,0.0)-change)#/f(xface0,0.0)
diff2 = (f(xface0,dt)-change2)#/f(xface0,dt)
diff3 = (f(xface0,2.0*dt)-change3)#/f(xface0,2.0*dt)
diff4 = (f(xface0,3.0*dt)-change4)#/f(xface0,3.0*dt)

diff5 = (f(xface0,4.0*dt)-change5)#/f(xface0,4.0*dt)
diff6 = (f(xface0,5.0*dt)-change6)#/f(xface0,5.0*dt)
diff7 = (f(xface0,6.0*dt)-change7)#/f(xface0,6.0*dt)
diff8 = (f(xface0,7.0*dt)-change8)#/f(xface0,7.0*dt)
diff9 = (f(xface0,8.0*dt)-change9)#/f(xface0,8.0*dt)


plt.xlim(0,20)
#plt.ylim(-0.5,0.5)
plt.xlabel("X [m]")
plt.ylabel("Difference [kg/m^3/s]")
plt.plot(x0,diff)
plt.plot(x0,diff2)
plt.plot(x0,diff3)
plt.plot(x0,diff4)

#------------------------------------------------------------------------------

plt.figure(2)

time = np.array([0.0, dt, 2.0*dt, 3.0*dt, 4.0*dt, 5.0*dt, 6.0*dt, 7.0*dt, 8.0*dt])
total_error = np.array([np.sum(abs(diff)),np.sum(abs(diff2)),np.sum(abs(diff3)),np.sum(abs(diff4)),np.sum(abs(diff5)),np.sum(abs(diff6)),np.sum(abs(diff7)),np.sum(abs(diff8)),np.sum(abs(diff9))])

plt.xlabel("Time [s]")
plt.ylabel("Total Absolute Difference in Density Flow [kg/m^3/s]")
plt.plot(time,total_error)


plt.show()


