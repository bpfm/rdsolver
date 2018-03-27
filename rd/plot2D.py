import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

n_points = 50

slide = 0

data = np.loadtxt("column.txt")

x = data[0:n_points,0]

rho0 = data[(slide+0)*n_points:(slide+1)*n_points,1]
rho1 = data[(slide+1)*n_points:(slide+2)*n_points,1]
rho2 = data[(slide+2)*n_points:(slide+3)*n_points,1]
rho3 = data[(slide+8)*n_points:(slide+9)*n_points,1]

pres0 = data[(slide+0)*n_points:(slide+1)*n_points,2]
pres1 = data[(slide+1)*n_points:(slide+2)*n_points,2]
pres2 = data[(slide+2)*n_points:(slide+3)*n_points,2]
pres3 = data[(slide+8)*n_points:(slide+9)*n_points,2]

x_vel0 = data[(slide+0)*n_points:(slide+1)*n_points,3]
x_vel1 = data[(slide+1)*n_points:(slide+2)*n_points,3]
x_vel2 = data[(slide+2)*n_points:(slide+3)*n_points,3]
x_vel3 = data[(slide+8)*n_points:(slide+9)*n_points,3]

plt.subplot(311)
plt.ylabel("Density [kg/m^3]")
#plt.plot(x,rho0)
#plt.plot(x,rho1)
#plt.plot(x,rho2)
plt.plot(x,rho3)
#plt.plot(x,rho1-rho0)
#plt.plot(x,rho2-rho1)

plt.subplot(312)
plt.ylabel("Pressure [N/m^2]")
#plt.plot(x,pres0)
#plt.plot(x,pres1)
#plt.plot(x,pres2)
plt.plot(x,pres3)
#plt.plot(x,pres1-pres0)
#plt.plot(x,pres2-pres1)

plt.subplot(313)
plt.xlabel("X [m]")
plt.ylabel("X Velocity [m/s]")
#plt.plot(x,x_vel0)
#plt.plot(x,x_vel1)
#plt.plot(x,x_vel2)
plt.plot(x,x_vel3)
#plt.plot(x,x_vel1-x_vel0)
#plt.plot(x,x_vel2-x_vel1)

plt.show()