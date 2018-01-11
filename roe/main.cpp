#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "constants.h"

#ifdef ONE_D
#include "centre1D.h"
#include "face1D.h"
#include "setup1D.cpp"
#include "io1D.cpp"
#endif

#ifdef THREE_D
#include "centre3D.h"
#include "face3D.h"
#include "setup3D.cpp"
#include "io3D.cpp"
#endif


using namespace std;

int main(){

        int i,j,k,l=0;                                          // ******* decalare varaibles and vectors ******
        double dx,dt,t=0.0,c_initial;                           // dx = space step,dt = timestep,t = time,cfl = cfl condition,c_initial = initial max
        double next_time=0.0;                                   // t_tot = total time,next_time = time of next snapshot
        centre new_centre;                                      // new_centre = temporary centre to be added to vector of vertices
        face new_x_face,new_y_face,new_z_face;                  // new_face = temporary face to be added to vector of faces
        vector<centre> points;                                  // points = vector of vertices
        vector<face> faces;                                     // faces = vector of faces
        vector<centre>::iterator it_vert;                       // it_vert = iterator for centre vector
        vector<face>::iterator it_face;                         // it_face = iterator for face vector
        double total_density,next_dt,possible_dt;               // total_density = total density in box
        int current_point=0;

        dx = SIDE_LENGTH/double(N_POINTS);                             // calculate cell width
        next_dt = t_tot;

        /****** Setup initial conditions of one dimensional tube ******/

        cout << "Building grid of cells ..." << endl;

#ifdef ONE_D
        for(i=0;i<N_POINTS;i++){
                new_centre = setup_centre(N_POINTS,i,dx);               // call centre setup routine
                points.push_back(new_centre);                           // add new centre to vector of all centres
                points[i].calc_next_dt(dx,cfl,possible_dt);             // check dt is min required by cfl
                if(possible_dt<next_dt){next_dt=possible_dt;}
        }
#endif

#ifdef THREE_D
        for(i=0;i<N_POINTS;i++){
                cout << "setting up x =\t" << (i+0.5)*dx << endl;
                for(j=0;j<N_POINTS;j++){
                        for(k=0;k<N_POINTS;k++){
                                new_centre = setup_centre(N_POINTS,i,j,k,dx);                   // call centre setup routine
                                points.push_back(new_centre);                                   // add new centre to vector of all centres
                                points[current_point].calc_next_dt(dx,cfl,possible_dt);         // check dt is min required by cfl
                                if(possible_dt<next_dt){next_dt=possible_dt;}
                                current_point += 1;
                        }
                }
        }
#endif

        /****** Setup system of faces ******/

        cout << "Assigning cells to faces ..." << endl;

#ifdef ONE_D
        for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
                new_x_face = setup_face(i, points);
                faces.push_back(new_x_face);                                      // add new face to vector of all faces
        }
#endif

#ifdef THREE_D
        for(it_vert=points.begin(),i=0;it_vert<points.end();it_vert++,i++){
                cout << "setting up faces at\t" << points[i].get_x() << "\t" << points[i].get_y() << "\t" << points[i].get_z() << endl;
                new_x_face = setup_face(i,dx,points,"x");
                new_y_face = setup_face(i,dx,points,"y");
                new_z_face = setup_face(i,dx,points,"z");
                faces.push_back(new_x_face);                                      // add new faces to vector of all faces
                faces.push_back(new_y_face);
                faces.push_back(new_z_face);
        }
#endif
        /****** Loop over time until total time t_tot is reached ******/

        ofstream density_map, pressure_map, velocity_map, density_slice, du_file;

        open_files(density_map, pressure_map, velocity_map, density_slice, du_file);           // open output files

        cout << "Evolving fluid ..." << endl;

        while(t<t_tot){

                dt = next_dt;                                                   // set timestep based oncaclulation from previous timestep

                dt = 0.001;

                total_density = 0.0;                                            // reset total density counter

                if(t>=next_time){                                               // write out densities at given interval
                        next_time=next_time+t_tot/float(N_SNAP);
                        if(next_time>t_tot){next_time=t_tot;}
                        output_state(density_map, pressure_map, velocity_map, density_slice, du_file, points, t, dt, dx);
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
                l+=1;                                                                           // increment step number
                cout << l << "\t" << t << endl;
        }

        output_state(density_map, pressure_map, velocity_map, density_slice, du_file, points, t, dt, dx);      // write out final state

        close_files(density_map, pressure_map, velocity_map, density_slice, du_file);

        return 0;
}
