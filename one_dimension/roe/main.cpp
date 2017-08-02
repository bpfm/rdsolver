#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "centre.h"
#include "face.h"
#include "setup.cpp"

using namespace std;

extern centre setup(int n_points, int i, float dx, centre new_centre, int ic);

int main(){

	int i,j;						// ******* decalare varaibles and vectors ******
	int n_points=40;					// n_points = number of vertices
	int centre_id_0,centre_id_1;				// centre_id_0 and centre_id_1 = number of
	double dx,dt,t=0.0,cfl,c_initial;			// dx = space step,dt = timestep,t = time,cfl = cfl condition,c_initial = initial max
	double t_tot=0.01,next_time=0.0;			// t_tot = total time,next_time = time of next snapshot
	centre new_centre;					// new_centre = temporary centre to be added to vector of vertices
	face new_face;						// new_face = temporary face to be added to vector of faces
	centre *centre_0,*centre_1,*centre_00,*centre_11;	// *centre_0 and *centre_1 = pointers to vertices
	vector<centre> points;					// points = vector of vertices
	vector<face> centers;					// centers = vector of faces
	vector<centre>::iterator it_vert;			// it_vert = iterator for centre vector
	vector<face>::iterator it_face;				// it_face = iterator for face vector
	double total_density,next_dt,possible_dt;		// total_density = total density in box

	dx = 40.0/double(n_points);				// calculate face width
	cfl = 0.1;						// set CFL condition
	next_dt = 1.0;

	cout << "Initial Timestep chosen as " << next_dt << " s" << endl;

	/****** Setup initial conditions of one dimensional tube ******/

	for (i=0;i<n_points;i++){
		new_centre = setup_centre(n_points,i,dx,new_centre,0);	// call centre setup routine (0 = shock tube, 1 = sine wave)
		points.push_back(new_centre);				// add new centre to vector of all vertices
		points[i].calc_next_dt(dx,cfl,possible_dt);		// check dt is min required by cfl
		if(possible_dt<next_dt){next_dt=possible_dt;}
	}

	/****** Setup system of faces ******/

	for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
		centre_id_0 = i % n_points;					// setup preiodic boundary
		centre_id_1 = (i+1) % n_points;

		//cout << i << " " << centre_id_0-1 << " " << centre_id_0 << " " << centre_id_1 << " " << centre_id_1+1 << endl;

		centre_0 = &points[centre_id_0];				// setup pointers to lower and upper centre
		centre_1 = &points[centre_id_1];
		centre_00 = &points[centre_id_0-1];				// add neighbouring vertices for flux limiter
		centre_11 = &points[centre_id_1+1];

		if(centre_id_0 == 0){centre_00 = &points[n_points-1];}
		if(centre_id_1 == n_points-1){centre_11 = &points[0];}

		new_face.set_centre_0(centre_0);				// pass these pointers to the face
		new_face.set_centre_1(centre_1);
		new_face.set_centre_00(centre_00);
		new_face.set_centre_11(centre_11);

		centers.push_back(new_face);					// add new face to vector of all faces
	}

	/****** Loop over time until total time t_tot is reached ******/

	ofstream density_map;						        // open output file
        ofstream pressure_map;

	density_map.open("density.txt");
        pressure_map.open("pressure.txt");

	while(t<t_tot){

		dt = next_dt;
		dt = 0.000001;
		total_density = 0.0;						// reset total density counter

		if(t>=next_time){						// write out densities at given interval
			next_time=next_time+(t_tot/10.0);
			for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
				density_map << points[i].get_x() << "\t" << points[i].get_mass_density() << endl;
                                pressure_map << points[i].get_x() << "\t" << points[i].get_pressure() << endl;
				total_density += points[i].get_mass_density()*dx;
			}
			cout << "*************************************" << endl;		// right out time and total density to terminal
			cout << "time " << t << " -> total mass = " << total_density  << " time step =  " << dt << endl;
			density_map << " " << endl;
                        pressure_map << " " << endl;
		}

		for(it_face=centers.begin(),i=0;it_face<centers.end();it_face++,i++){		// loop over all faces
			centers[i].construct_state(dx,dt,t);					// calculate flux through boundary
		}

		for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){		// loop over all vertices
			points[i].update_u_variables();						// update the u variables with the collected du
			points[i].prim_to_con();						// convert these to their corresponding conserved
                        points[i].recalculate_pressure();                                       // caclulate pressure from new conserved values
			points[i].con_to_prim();						// convert back to guarentee correct values are used
			points[i].setup_f_variables();						// set flux variables with new values
			points[i].reset_du();							// reset du value to zero for next timestep
			points[i].calc_next_dt(dx,cfl,possible_dt);				// calculate next timestep
			if(possible_dt<next_dt){next_dt=possible_dt;}
		}
		t+=dt;																	// increment time
	}
	density_map.close();
        pressure_map.close();
	return 0;
}
