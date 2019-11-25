import numpy as np
import matplotlib.pyplot as plt

xlow, ylow, xhigh, yhigh = np.loadtxt("output.tri",skiprows=0,max_rows=1)
nx, ny = np.loadtxt("output.tri",skiprows=1,max_rows=1)
ndata = 2*np.loadtxt("output.tri",skiprows=2,max_rows=1)

npoint = int(ndata/2)

data = np.loadtxt("output.tri",skiprows=3,max_rows=int(ndata))

point  = np.zeros(shape=(int(ndata/2),2))
offset = np.zeros(shape=(int(ndata/2),2))

x = np.zeros(shape=npoint)
y = np.zeros(shape=npoint)

for i in range(0,int(ndata),2):
	point[int(i/2),:] = data[i]

for i in range(1,int(ndata),2):
	offset[int((i-1)/2),:] = data[i]

for j in range(0,npoint):
	x[j] = point[j,0] + offset[j,0]*(xhigh - xlow)
	y[j] = point[j,1] + offset[j,1]*(yhigh - ylow)

# print(point)

ntri = np.loadtxt("output.tri",skiprows=3+int(ndata),max_rows=1)

data_tri = np.loadtxt("output.tri",skiprows=5+int(ndata),max_rows=int(ntri))

# print(data_tri)

v0 = np.zeros(shape=len(data_tri))
v1 = np.zeros(shape=len(data_tri))
v2 = np.zeros(shape=len(data_tri))

def connectpoints(x,y,p1,p2):
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1], y[p2]
    if (x1 > 2.0 and x1 < 4.0) and (y1 > 2.0 and y1 < 4.0) and (x2 > 2.0 and x2 < 4.0) and (y2 > 2.0 and y2 < 4.0):
    	plt.plot([x1,x2],[y1,y2],'k-')
    elif (x1 > 2.0 and x1 < 4.0) and (y1 > 2.0 and y1 < 4.0) or (x2 > 2.0 and x2 < 4.0) and (y2 > 2.0 and y2 < 4.0):
    	plt.plot([x1,x2],[y1,y2],'g-')


for j in range(int(ntri)):
	v0 = int(data_tri[j][0])
	v1 = int(data_tri[j][1])
	v2 = int(data_tri[j][2])

	connectpoints(x,y,v0,v1)
	connectpoints(x,y,v0,v2)
	connectpoints(x,y,v1,v2)

plt.scatter(x,y)
plt.xlim(1.5,4.5)
plt.ylim(1.5,4.5)

plt.vlines(2.0,ylow - ny,yhigh*(ny+1),color='red')
plt.vlines(4.0,ylow - ny,yhigh*(ny+1),color='red')
plt.hlines(2.0,xlow - nx,xhigh*(nx+1),color='red')
plt.hlines(4.0,xlow - nx,xhigh*(nx+1),color='red')


plt.show()