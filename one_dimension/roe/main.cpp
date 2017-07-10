#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "vertex.h"
#include "cell.h"

using namespace std;

int main(){

	int i,j;							// ******* decalare varaibles and vectors ******
	int n_points=100;					// n_points = number of vertices
	int vertex_id_0,vertex_id_1;		// vertex_id_0 and vertex_id_1 = number of
	double dx,dt,t=0.0,cfl,c_initial;	// dx = space step,dt = timestep,t = time,cfl = cfl condition,c_initial = initial max sound speed
	double t_tot=0.001,next_time=0.0;
	vertex new_vertex;					// new_vertex = temporary vertex to be added to vector of vertices
	cell new_cell;						// new_cell = temporary cell to be added to vector of cells
	vertex *vertex_0,*vertex_1;			// *vertex_0 and *vertex_1 = pointers to vertices
	vector<vertex> points;				// points = vector of vertices
	vector<cell> centers;				// centers = vector of cells
	vector<vertex>::iterator it_vert;	// it_vert = iterator for vertex vector
	vector<cell>::iterator it_cell;		// it_cell = iterator for cell vector
	float total_density;

	dx = 20.0/double(n_points);

	cfl = 0.001;
    c_initial = 374.17;			// units m/s
    dt = cfl*(dx/c_initial);

    //srand(static_cast <unsigned> (time(0)));

	/****** Setup initial conditions of one dimensional tube ******/

	for (i=0;i<n_points;i++){
			new_vertex.set_x(float(i)*dx);

			if(i<n_points/2){
				new_vertex.set_mass_density(1.0);						// units kg/m^3
				new_vertex.set_velocity(0.0);							// units m/s
				new_vertex.set_pressure(100000.0);						// units N/m^2
			}else{
				new_vertex.set_mass_density(0.125);						// units kg/m^3
				new_vertex.set_velocity(0.0);							// units m/s
				new_vertex.set_pressure(1000.0);						// units N/m^2
			}
			new_vertex.setup_energy_density();
			new_vertex.con_to_prim();
			new_vertex.setup_f_variables();
			new_vertex.reset_du();
			points.push_back(new_vertex);			// add new vertex to vector of all vertices
		}

	/****** Setup system of cells ******/

	for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
		vertex_id_0 = i % n_points;				// setup preiodic boundary
		vertex_id_1 = (i+1) % n_points;

		vertex_0 = &points[vertex_id_0];		// setup pointers to lower and upper vertex
		vertex_1 = &points[vertex_id_1];

		new_cell.set_vertex_0(vertex_0);		// pass these pointers to the cell
		new_cell.set_vertex_1(vertex_1);

		centers.push_back(new_cell);			// add new cell to vector of all cells
	}

	/****** Loop over total time t_tot ******/

	ofstream density_map;
	density_map.open("density.txt");

	while(t<t_tot){

		total_density = 0.0;		// reset total density counter

		if(t>=next_time){			// write out densities at given interval

			next_time=next_time+(t_tot/10.0);

			for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
				density_map << points[i].get_x() << "\t" << points[i].get_mass_density() << endl;
				total_density += points[i].get_mass_density()*dx;
			}

			cout << "time " << t << " -> total mass = " << total_density << endl;
			density_map << " " << endl;

		}

		for(it_cell=centers.begin(),i=0;it_cell<centers.end();it_cell++,i++){	// loop over all cells
			centers[i].construct_state(dx,dt);									// calculate flux through boundary
		}

		for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){		// loop over all vertices
			//cout << points[i].get_u0() << " " << points[i].get_u1() << " " << points[i].get_u2() << endl;
			points[i].update_u_variables();										// update the u variables with the collected du
			points[i].prim_to_con();											// convert these to their corresponding conserved
			//points[i].con_to_prim();
			points[i].setup_f_variables();										// set flux variables with new values
			points[i].reset_du();												// reset du value to zero for next timestep
		}

		t+=dt;								// increment time

	}

	density_map.close();
	return 0;
}