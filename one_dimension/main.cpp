#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "vertex.h"
#include "cell.h"


using namespace std;

int main(){

	int i,j;
	int n_points=10,n_steps=10;
	int vertex_id_0,vertex_id_1;
	float dx,dt,t_total=1.0,t=0.0;
	vertex *vertex_0,*vertex_1;
	vertex new_vertex;
	cell new_cell;
	vector<vertex> points;
	vector<cell> centers;
	vector<vertex>::iterator it_vert;
	vector<cell>::iterator it_cell;

	dt=t_total/float(n_steps);
	dx=10.0/float(n_points);

	/****** Setup initial conditions of one dimensional tube ******/

	for (i=0;i<n_points;i++){
		new_vertex.set_x(float(i));
		if(i<n_points/2){
			new_vertex.set_mass_density(1.0+0.0000001*float(i));	// units kg/m^3
			new_vertex.set_velocity(0.0000001*float(i));			// units m/s
			new_vertex.set_pressure(100);							// units kN/m^2
		}else{
			new_vertex.set_mass_density(0.125+0.0000001*float(i));	// units kg/m^3
			new_vertex.set_velocity(-0.0000001*float(i));			// units m/s
			new_vertex.set_pressure(10);							// units kN/m^2
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

		//cout << "vertex_id " << vertex_id_0 << " " << vertex_id_1 << endl;

		vertex_0 = &points[vertex_id_0];		// setup pointers to lower and upper vertex
		vertex_1 = &points[vertex_id_1];

		new_cell.set_vertex_0(vertex_0);		// pass these pointers to the cell
		new_cell.set_vertex_1(vertex_1);

		centers.push_back(new_cell);			// add new cell to vector of al cells
	}

	/****** Loop over n_steps dt ******/

	for(j=0;j<n_steps;j++){
		
		for(it_cell=centers.begin(),i=0;it_cell<centers.end();it_cell++,i++){
			centers[i].distribute_residual(dx,dt);
		}
		for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
			points[i].update_u_variables();
			points[i].prim_to_con();
			points[i].reset_du();
			points[i].con_to_prim();
			cout << dt << " " << points[i]
			.get_mass_density() << endl;
		}
		t+=dt;
	}
	return 0;
}