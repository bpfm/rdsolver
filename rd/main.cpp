#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include "constants.h"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#include "setup2D.cpp"
#include "io2D.cpp"
#endif

using namespace std;

int main(){

        int i,j,k,l=0;                                          // ******* decalare varaibles and vectors ******
        double DX,DY,DT,T=0.0;                                     // DX = space step,DT = timestep,t = time,CFL = CFL condition
        double NEXT_TIME=0.0;                                   // T_TOT = total time,NEXT_TIME = time of next snapshot
        VERTEX NEW_VERTEX;
        TRIANGLE NEW_TRIANGLE;
        vector<VERTEX> X_POINTS;
        vector<vector<VERTEX> > POINTS;
        vector<TRIANGLE> X_MESH;
        vector<vector<TRIANGLE> > MESH;
        vector<VERTEX>::iterator IT_VERT;
        vector<TRIANGLE>::iterator IT_TRIANGLE;
        double TOTAL_DENSITY,NEXT_DT,POSSIBLE_DT;               // TOTAL_DENSITY = total density in box

        NEXT_DT = T_TOT;

        /****** Setup initial conditions of one dimensional tube ******/

        cout << "Building grid of vertices ..." << endl;

#ifdef TWO_D
        for(j=0;j<N_POINTS;j++){
                for(i=0;i<N_POINTS;i++){
                        NEW_VERTEX = setup_vertex(N_POINTS,i,j,DX,DY);               // call VERTEX setup routine
                        X_POINTS.push_back(NEW_VERTEX);                           // add new VERTEX to vector of vertices in this row
                        //POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);                      // check dt is min required by CFL
                        //if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                }
                POINTS.push_back(X_POINTS);
                X_POINTS.clear();
        }
#endif

        /****** Setup system of MESH ******/

         cout << "Assigning vertices to triangles ..." << endl;

        for(j=0;j<2*(N_POINTS);j++){
                for(i=0;i<N_POINTS;i++){
                        NEW_TRIANGLE = setup_triangle(i,j,POINTS);
                        X_MESH.push_back(NEW_TRIANGLE);
                }
                MESH.push_back(X_MESH);
                X_MESH.clear();
        }

        cout << "Finished triangle setup" << endl;

        /****** Loop over time until total time T_TOT is reached ******/

        ofstream DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE;

        open_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);           // open output files

        cout << "Evolving fluid ..." << endl;

        cout << "Mesh Size =\t" << MESH[0].size() << '\t' << MESH.size() << endl;

        while(T<T_TOT){

                //cout << "STEP =\t" << l << endl;

                DT = NEXT_DT;                                                   // set timestep based oncaclulation from previous timestep

                DT = 0.000001;

                TOTAL_DENSITY = 0.0;                                            // reset total density counter

                cout << "time =\t" << T << endl;

                if(T >= 0.9999999*NEXT_TIME){                                             // write out densities at given interval
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DY);
                }

                for(j=0;j<2*(N_POINTS);j++){
                        for(i=0;i<N_POINTS;i++){                               // loop over all MESH
                                MESH[j][i].calculate_change(T, DT);          // calculate flux through TRIANGLE
                        }
                }

                NEXT_DT = T_TOT - (T + DT);        // set next timestep to max possible value (time remaining to end)

                for(j=0;j<N_POINTS;j++){
                        for(i=0;i<N_POINTS;i++){             // loop over all vertices
                                POINTS[j][i].update_u_variables();                                         // update the u variables with the collected du
                                POINTS[j][i].con_to_prim();                                                // convert these to their corresponding conserved
                                POINTS[j][i].recalculate_pressure();                                       // caclulate pressure from new conserved values
                                POINTS[j][i].prim_to_con();                                                // convert back to guarentee correct values are used
                                POINTS[j][i].reset_du();                                                   // reset du value to zero for next timestep
                                //POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);                                // calculate next timestep
                                //if(POSSIBLE_DT<NEXT_DT){NEXT_DT = POSSIBLE_DT;}
                        }
                }
                T+=DT;                                                                          // increment time
                l+=1;                                                                           // increment step number

        }

         output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DX);      // write out final state

         close_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);

        return 0;
}
