#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "vertex.h"
#include "cell.h"
#include "setup.cpp"

using namespace std;

extern vertex setup(int n_points, int i, float dx, vertex new_vertex, int ic);

int main(){

	int i,j;						// ******* decalare varaibles and vectors ******
	int n_points=200;					// n_points = number of vertices
	int vertex_id_0,vertex_id_1;				// vertex_id_0 and vertex_id_1 = number of
	double dx,dt,t=0.0,cfl,c_initial;			// dx = space step,dt = timestep,t = time,cfl = cfl condition,c_initial = initial max
	double t_tot=0.001,next_time=0.0;			// t_tot = total time,next_time = time of next snapshot
	vertex new_vertex;					// new_vertex = temporary vertex to be added to vector of vertices
	cell new_cell;						// new_cell = temporary cell to be added to vector of cells
	vertex *vertex_0,*vertex_1,*vertex_00,*vertex_11;	// *vertex_0 and *vertex_1 = pointers to vertices
	vector<vertex> points;					// points = vector of vertices
	vector<cell> centers;					// centers = vector of cells
	vector<vertex>::iterator it_vert;			// it_vert = iterator for vertex vector
	vector<cell>::iterator it_cell;				// it_cell = iterator for cell vector
	double total_density,next_dt,possible_dt;		// total_density = total density in box

	dx = 20.0/double(n_points);				// calculate cell width
	cfl = 0.1;						// set CFL condition
	next_dt = 1.0;

	cout << "Initial Timestep chosen as " << next_dt << " s" << endl;

	/****** Setup initial conditions of one dimensional tube ******/

	for (i=0;i<n_points;i++){
		new_vertex = setup_vertex(n_points,i,dx,new_vertex,0);	// call vertex setup routine (0 = shock tube, 1 = sine wave)
		points.push_back(new_vertex);				// add new vertex to vector of all vertices
		points[i].calc_next_dt(dx,cfl,possible_dt);		// check dt is min required by cfl
		if(possible_dt<next_dt){next_dt=possible_dt;}
	}

	/****** Setup system of cells ******/

	for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
		vertex_id_0 = i % n_points;					// setup preiodic boundary
		vertex_id_1 = (i+1) % n_points;

		//cout << i << " " << vertex_id_0-1 << " " << vertex_id_0 << " " << vertex_id_1 << " " << vertex_id_1+1 << endl;

		vertex_0 = &points[vertex_id_0];				// setup pointers to lower and upper vertex
		vertex_1 = &points[vertex_id_1];
		vertex_00 = &points[vertex_id_0-1];				// add neighbouring vertices for flux limiter
		vertex_11 = &points[vertex_id_1+1];

		if(vertex_id_0 == 0){vertex_00 = &points[n_points-1];}
		if(vertex_id_1 == n_points-1){vertex_11 = &points[0];}

		new_cell.set_vertex_0(vertex_0);				// pass these pointers to the cell
		new_cell.set_vertex_1(vertex_1);
		new_cell.set_vertex_00(vertex_00);
		new_cell.set_vertex_11(vertex_11);

		centers.push_back(new_cell);					// add new cell to vector of all cells
	}

	/****** Loop over time until total time t_tot is reached ******/

	ofstream density_map;							// open output file
	density_map.open("density.txt");					

	while(t<t_tot){

		dt = next_dt;
		total_density = 0.0;						// reset total density counter

		if(t>=next_time){						// write out densities at given interval
			next_time=next_time+(t_tot/10.0);
			for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
				density_map << points[i].get_x() << "\t" << points[i].get_mass_density() << endl;
				total_density += points[i].get_mass_density()*dx;
			}
			cout << "*************************************" << endl;		// right out time and total density to terminal
			cout << "time " << t << " -> total mass = " << total_density  << " time step =  " << dt << endl;
			density_map << " " << endl;
		}

		for(it_cell=centers.begin(),i=0;it_cell<centers.end();it_cell++,i++){		// loop over all cells
			centers[i].construct_state(dx,dt,t);					// calculate flux through boundary
		}

		for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){		// loop over all vertices
			points[i].update_u_variables();						// update the u variables with the collected du
			points[i].prim_to_con();						// convert these to their corresponding conserved
			points[i].con_to_prim();						// convert back to guarentee correct value are used
			points[i].setup_f_variables();						// set flux variables with new values
			points[i].reset_du();							// reset du value to zero for next timestep
			points[i].calc_next_dt(dx,cfl,possible_dt);				// calculate next timestep
			if(possible_dt<next_dt){next_dt=possible_dt;}
		}
		t+=dt;																	// increment time
	}
	density_map.close();
	return 0;
}
