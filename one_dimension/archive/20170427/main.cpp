#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "cell.h"
#include "euler.cpp"
#include "integrate.cpp"
#include "residual.cpp"

using namespace std;

extern float approximation(float x,float lower,float upper, float lower_value,float upper_value);

int main(){

	int i,j,
	n_cell=10,n_points=10;
	float new_x,new_lower,new_upper,new_mass,new_momentum,new_energy;
	vector<cell> domain,old_domain;
	vector<cell>::iterator it1;
	cell new_cell;
	float Q_h[n_points],x_pos;

	for(i=0;i<n_cell;i++){

		new_x = float(i);
		new_lower = new_x - 0.5;
		new_upper = new_x + 0.5;
		new_mass = 1.0;
		new_momentum = 1.0;
		new_energy = 1.0;

		new_cell.set_x(new_x);
		new_cell.set_vertex_lower(new_lower);
		new_cell.set_vertex_upper(new_upper);
		new_cell.set_mass(new_mass);
		new_cell.set_momentum(new_momentum);
		new_cell.set_energy(new_energy);
		new_cell.con_to_prim();
		new_cell.setup_u_variables();
		new_cell.setup_fluxes();

		domain.push_back(new_cell);
	}

	for(it1=domain.begin(),i=0;it1<domain.end();it1++,i++){
		for(j=0;j<n_points;j++){
			x_pos = domain[i].get_vertex_lower()+(domain[i].get_vertex_upper()-domain[i].get_vertex_lower())*float(j)/float(n_points);
			Q_h[j] = approximation(x_pos,domain[i].get_vertex_lower(),domain[i].get_vertex_upper(),domain[i].get_u0(),domain[i].get_u0());
		}
	}
	
}