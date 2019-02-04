import matplotlib.pyplot as plt
import numpy as np

#x,density = np.loadtxt("density.txt", usecols = (0,1))

n_snap = 5

n_side = 50
side_size = 50.0
n_points = n_side*n_side*n_side
dx = side_size/float(n_side)

data = np.loadtxt("density.txt")
x0,y0,z0,rho0 = data[0:n_points,0]-0.50*dx, data[0:n_points,1]-0.50*dx, data[0:n_points,2]-0.50*dx, data[0:n_points,3]
x1,y1,z1,rho1 = data[n_points:2*n_points,0]-0.50*dx, data[n_points:2*n_points,1]-0.50*dx, data[n_points:2*n_points,2]-0.50*dx, data[n_points:2*n_points,3]
x2,y2,z2,rho2 = data[2*n_points:3*n_points,0]-0.50*dx, data[2*n_points:3*n_points,1]-0.50*dx, data[2*n_points:3*n_points,2]-0.50*dx, data[2*n_points:3*n_points,3]
x3,y3,z3,rho3 = data[3*n_points:4*n_points,0]-0.50*dx, data[3*n_points:4*n_points,1]-0.50*dx, data[3*n_points:4*n_points,2]-0.50*dx, data[3*n_points:4*n_points,3]
x4,y4,z4,rho4 = data[4*n_points:5*n_points,0]-0.50*dx, data[4*n_points:5*n_points,1]-0.50*dx, data[4*n_points:5*n_points,2]-0.50*dx, data[4*n_points:5*n_points,3]

x_range = np.zeros(shape = n_side)
y_range = np.zeros(shape = n_side)
z_range = np.zeros(shape = n_side)

for j in range(n_snap):
        density_grid = np.zeros(shape = (n_side,n_side,n_side))
        for i in range(n_points):
                #print(x0[i])
                x_bin = int(x0[i]*n_side/side_size)
                y_bin = int(y0[i]*n_side/side_size)
                z_bin = int(z0[i]*n_side/side_size)
                if j==0:
                        density_grid[x_bin,y_bin,z_bin] = rho0[i]
                if j==1:
                        density_grid[x_bin,y_bin,z_bin] = rho1[i]
                if j==2:
                        density_grid[x_bin,y_bin,z_bin] = rho2[i]
                if j==3:
                        density_grid[x_bin,y_bin,z_bin] = rho3[i]
                if j==4:
                        density_grid[x_bin,y_bin,z_bin] = rho4[i]
                        
                x_range[x_bin] = x0[i]
                y_range[y_bin] = y0[i]
                z_range[z_bin] = z0[i]

        density_collapse_z = np.sum(density_grid,1)/float(n_side)

        density_collapse_yz = np.sum(density_collapse_z,1)/float(n_side)

        #plt.figure(100+j)

        #plt.pcolor(x_range,y_range,density_collapse_z,cmap=plt.cm.Oranges)

        #plt.xlabel("x [m]")
        #plt.ylabel("y [m]")
        #plt.xlim([20,30])
        #plt.ylim([0.001,1.5])

        plt.figure(110)

        plt.plot(x_range,density_collapse_yz)

        plt.xlabel("x [m]")
        plt.ylabel("density [kg/m^3]")

plt.show()
