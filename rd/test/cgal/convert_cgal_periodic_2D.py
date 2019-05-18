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



ntri = np.loadtxt("output.tri",skiprows=3+int(ndata),max_rows=1)

data_tri = np.loadtxt("output.tri",skiprows=5+int(ndata),max_rows=int(ntri))

v0 = np.zeros(shape=len(data_tri))
v1 = np.zeros(shape=len(data_tri))
v2 = np.zeros(shape=len(data_tri))

main_xlow  = xlow + xhigh
main_xhigh = xlow + 2.0*xhigh

main_ylow  = ylow + yhigh
main_yhigh = ylow + 2.0*yhigh

f = open('triangles.txt','w')

# f.write(str(int(ntri)) + "\n") 

for j in range(int(ntri)):
	v0 = int(data_tri[j][0])
	v1 = int(data_tri[j][1])
	v2 = int(data_tri[j][2])
	x0 = x[v0]
	x1 = x[v1]
	x2 = x[v2]
	y0 = y[v0]
	y1 = y[v1]
	y2 = y[v2]

	# print(main_xlow,main_xhigh,main_ylow,main_yhigh)

	# check if any vertex of triangle is within main region
	if ((x0 > main_xlow and x0 < main_xhigh) and (y0 > main_ylow and y0 < main_yhigh))\
	or ((x1 > main_xlow and x1 < main_xhigh) and (y1 > main_ylow and y1 < main_yhigh))\
	or ((x2 > main_xlow and x2 < main_xhigh) and (y2 > main_ylow and y2 < main_yhigh)):
		f.write(str(v0) + "\t" + str(v1) + "\t" + str(v2) + "\n")

