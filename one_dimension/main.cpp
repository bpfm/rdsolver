#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "vertex.h"
#include "cell.h"


using namespace std;

int main(){

	int i;
	int n_points=10;
	int vertex_id_0,vertex_id_1;
	vertex *vertex_0,*vertex_1;
	vertex new_vertex;
	cell new_cell;
	vector<vertex> points;
	vector<cell> centers;
	vector<vertex>::iterator it_vert;
	vector<cell>::iterator it_cell;

	for (i=0;i<n_points;i++){
		new_vertex.set_x(float(i));
		if(i<n_points/2){
			new_vertex.set_mass_density(1.0);	// units kg/m^3
			new_vertex.set_velocity(1.0);	// units m/s
			new_vertex.set_pressure(100);	// units kN/m^2
		}else{
			new_vertex.set_mass_density(0.125);	// units kg/m^3
			new_vertex.set_velocity(-1.0);	// units m/s
			new_vertex.set_pressure(10);	// units kN/m^2
		}

		new_vertex.setup_energy_density();
		new_vertex.setup_u_variables();
		new_vertex.setup_f_variables();

		points.push_back(new_vertex);
		cout << new_vertex.get_x() << " " << new_vertex.get_mass_density() << endl;
	}

	for (it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
		vertex_id_0 = i % n_points;			// setup preiodic boundary
		vertex_id_1 = i+1 % n_points;

		vertex_0 = &points[vertex_id_0];	// setup pointers to lower and upper vertex
		vertex_1 = &points[vertex_id_1];

		new_cell.set_vertex_0(vertex_0);	// pass these pointers to the cell
		new_cell.set_vertex_1(vertex_1);

		centers.push_back(new_cell);

		cout << new_cell.get_element_residual() << endl;

	}

	return 0;
}