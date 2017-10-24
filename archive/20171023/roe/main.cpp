#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "constants.h"
#include "centre.h"
#include "face.h"
#include "io.cpp"
#include "setup.cpp"

using namespace std;

int main(){

        int i,j;                                                // ******* decalare varaibles and vectors ******
        double dx,dt,t=0.0,c_initial;                           // dx = space step,dt = timestep,t = time,cfl = cfl condition,c_initial = initial max
        double next_time=0.0;                                   // t_tot = total time,next_time = time of next snapshot
        centre new_centre;                                      // new_centre = temporary centre to be added to vector of vertices
        face new_face;                                          // new_face = temporary face to be added to vector of faces
        vector<centre> points;                                  // points = vector of vertices
        vector<face> faces;                                     // faces = vector of faces
        vector<centre>::iterator it_vert;                       // it_vert = iterator for centre vector
        vector<face>::iterator it_face;                         // it_face = iterator for face vector
        double total_density,next_dt,possible_dt;               // total_density = total density in box

        dx = 50.0/double(n_points);                             // calculate face width

        /****** Setup initial conditions of one dimensional tube ******/

        for (i=0;i<n_points;i++){
                new_centre = setup_centre(n_points,i,dx);    // call centre setup routine (0 = shock tube, 1 = sine wave)
                points.push_back(new_centre);                           // add new centre to vector of all vertices
                points[i].calc_next_dt(dx,cfl,possible_dt);             // check dt is min required by cfl
                if(possible_dt<next_dt){next_dt=possible_dt;}
        }

        /****** Setup system of faces ******/

        for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
                new_face = setup_face(i, points);
                faces.push_back(new_face);                                      // add new face to vector of all faces
        }

        /****** Loop over time until total time t_tot is reached ******/

        ofstream density_map, pressure_map, velocity_map, du_file;

        open_files(density_map, pressure_map, velocity_map, du_file);           // open output files

        j=1;

        while(t<t_tot){

                dt = next_dt;                                                   // set timestep based oncaclulation from previous timestep

                total_density = 0.0;                                            // reset total density counter

                if(t>=next_time){                                               // write out densities at given interval
                        next_time=next_time+0.1*t_tot;
                        if(next_time>t_tot){next_time=t_tot;}
                        output_state(density_map, pressure_map, velocity_map, du_file, points, t, dt, dx);
                }

                for(it_face=faces.begin(),i=0;it_face<faces.end();it_face++,i++){               // loop over all faces
                        faces[i].calculate_flux(dx,dt,t,du_file);                               // calculate flux through face
                }

                next_dt = t_tot - (t + dt);     // set next timestep to max possible value (time remaining to end)

                for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){             // loop over all vertices

                        points[i].update_u_variables();                                         // update the u variables with the collected du
                        points[i].con_to_prim();                                                // convert these to their corresponding conserved
                        points[i].recalculate_pressure();                                       // caclulate pressure from new conserved values
                        points[i].prim_to_con();                                                // convert back to guarentee correct values are used
                        points[i].setup_f_variables();                                          // set flux variables with new values
                        points[i].reset_du();                                                   // reset du value to zero for next timestep
                        points[i].calc_next_dt(dx,cfl,possible_dt);                             // calculate next timestep
                        if(possible_dt<next_dt){next_dt = possible_dt;}
                }

                t+=dt;                                                                          // increment time
                j+=1;                                                                           // increment step number
        }

        output_state(density_map, pressure_map, velocity_map, du_file, points, t, dt, dx);      // write out final state

        close_files(density_map, pressure_map, velocity_map, du_file);

        return 0;
}
