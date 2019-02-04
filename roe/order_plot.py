import matplotlib.pyplot as plt
import numpy as np

N = np.array([50000,5000,500,200,100,50,20,10])

dx = np.array([0.001,0.01,0.1,0.25,0.5,1.0,2.5,5.0])
error_dx = np.array([0.0340,0.342,3.347,8.060,15.249,27.698,68.656,106.067])#/(dx*N)

plt.subplot(211)

plt.xlabel("dx [m]")
plt.ylabel("Error [kg]")

print(error_dx)

plt.loglog(dx,error_dx)
#plt.plot(dx,error_dx)

plt.subplot(212)

dt = np.array([0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05])#,0.1])
error_dt = np.array([15.249,15.249,15.248,15.243,15.224,15.124,14.998,13.967])#/(0.5*100)#,12.622])/0.5

plt.xlabel("dt [s]")
plt.ylabel("Error [kg]")

print(error_dx)

#plt.semilogx(dt,error_dt)
plt.plot(dt,error_dt)

plt.show()