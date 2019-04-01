import numpy as np

NX = 1000
NY = 100

NTOT = NX*NY

XMIN = -0.5
XMAX = 0.5

YMIN = -0.5
YMAX = 0.5

f = open('points.txt','w')

f.write(str(2) + ' not_rbox ' + str(NTOT) + ' D2' +'\n')
f.write(str(NTOT) + '\n')

for i in range(NX):
	X = float(i)/float(NX)*(XMAX - XMIN) + XMIN
	for j in range(NY):
		Y = float(j)/float(NY)*(YMAX - YMIN) + YMIN
		f.write(str(X)+"\t"+str(Y)+'\n')

f.close()
