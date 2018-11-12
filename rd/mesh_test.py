import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("mesh_pos.txt")

x0,y0,x1,y1,x2,y2,xcen,ycen = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]

plt.plot(x0,y0,'ro')
for i in range(0,100):
  x01 = (x0[i],x1[i])
  y01 = (y0[i],y1[i])
  x02 = (x0[i],x2[i])
  y02 = (y0[i],y2[i])
  x12 = (x1[i],x2[i])
  y12 = (y1[i],y2[i])
  plt.plot(x01,y01,'b-')
  plt.plot(x02,y02,'r-')
  plt.plot(x12,y12,'g-')

plt.show()