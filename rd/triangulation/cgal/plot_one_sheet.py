import numpy as np
import matplotlib.pyplot as plt

xlow, ylow, xhigh, yhigh = np.loadtxt("output.tri",skiprows=0,max_rows=1)
nx, ny = np.loadtxt("output.tri",skiprows=1,max_rows=1)
ndata = np.loadtxt("output.tri",skiprows=2,max_rows=1)

data = np.loadtxt("output.tri",skiprows=3,max_rows=int(ndata))

x = np.zeros(shape=int(ndata))
y = np.zeros(shape=int(ndata))

for j in range(0,int(ndata)):
	x[j] = data[j,0]
	y[j] = data[j,1]

print(data)

ntri = np.loadtxt("output.tri",skiprows=3+int(ndata),max_rows=1)

data_tri = np.loadtxt("output.tri",skiprows=5+int(ndata),max_rows=int(ntri))

print(data_tri)

v0 = np.zeros(shape=len(data_tri))
v1 = np.zeros(shape=len(data_tri))
v2 = np.zeros(shape=len(data_tri))

def connectpoints(x,y,p1,p2):
    x1, x2 = x[p1], x[p2]
    y1, y2 = y[p1], y[p2]
    if (np.abs(x1 - x2) < 0.5*(xhigh - xlow)) and (np.abs(y1 - y2) < 0.5*(yhigh- ylow)):
    	plt.plot([x1,x2],[y1,y2],'k-')
    else:
    	plt.plot([x1,x2],[y1,y2],'g-')


for j in range(int(ntri)):
	v0 = int(data_tri[j][0])
	v1 = int(data_tri[j][1])
	v2 = int(data_tri[j][2])

	connectpoints(x,y,v0,v1)
	connectpoints(x,y,v0,v2)
	connectpoints(x,y,v1,v2)

# plt.scatter(data[:,0],data[:,1])
plt.show()