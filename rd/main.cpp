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

        int i, j, k, l = 0;                                          // ******* decalare varaibles and vectors ******
        double DX, DY, DT, T = 0.0;                                // DX = space step,DT = timestep,t = time,CFL = CFL condition
        double NEXT_TIME  = 0.0;                                   // T_TOT = total time,NEXT_TIME = time of next snapshot
        VERTEX                     NEW_VERTEX;
        TRIANGLE                   NEW_TRIANGLE;
        vector<VERTEX>             X_POINTS;
        vector<vector<VERTEX> >    POINTS;
        vector<TRIANGLE>           X_MESH;
        vector<vector<TRIANGLE> >  MESH;
        vector<VERTEX>::iterator   IT_VERT;
        vector<TRIANGLE>::iterator IT_TRIANGLE;
        double TOTAL_DENSITY, NEXT_DT, POSSIBLE_DT;               // TOTAL_DENSITY = total density in box

        NEXT_DT = T_TOT;

        /****** Setup initial conditions of one dimensional tube ******/

        cout << "Building grid of vertices" << endl;

#ifdef TWO_D
        for(j=0;j<N_POINTS;j++){
                for(i=0;i<N_POINTS;i++){
                        NEW_VERTEX = setup_vertex(N_POINTS,i,j,DX,DY);               // call VERTEX setup routine
                        X_POINTS.push_back(NEW_VERTEX);                              // add new VERTEX to vector of vertices in this row
                        //POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);                // check dt is min required by CFL
                        //if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                }
                POINTS.push_back(X_POINTS);
                X_POINTS.clear();
        }
#endif

        /****** Setup system of MESH ******/

         cout << "Assigning vertices to triangles" << endl;

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

        open_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);               // open output files

        cout << "Evolving fluid ..." << endl;
        cout << "Mesh Size =\t" << MESH[0].size() << '\t' << MESH.size() << endl;

        while(T<T_TOT){

                cout << "STEP =\t" << l << "\tTIME =\t" << T << endl;

                //DT = NEXT_DT;                                                     // set timestep based oncaclulation from previous timestep
                DT = 0.000001;

                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DY);
                }

                NEXT_DT = T_TOT - (T + DT);        // set next timestep to max possible value (time remaining to end)

                for(i=0;i<N_POINTS;i++){                                        // loop over all triangles in MESH
                        for(j=0;j<N_POINTS;j++){ 
                                MESH[i][j].calculate_first_half(T, DT);             // calculate flux through TRIANGLE
                        }
                }

                for(i=0;i<N_POINTS;i++){
                        for(j=0;j<N_POINTS;j++){                                     // loop over all vertices
                                POINTS[i][j].update_u_half();                        // update the half time state
                                POINTS[i][j].con_to_prim_half();
                                POINTS[i][j].recalculate_pressure_half();            // caclulate pressure from new conserved values
                                POINTS[i][j].reset_du_half();                             // reset du value to zero for next timestep
                        }
                }

                for(i=0;i<2*(N_POINTS);i++){                                         // loop over all triangles in MESH
                        for(j=0;j<N_POINTS;j++){
                                MESH[i][j].calculate_second_half(T, DT);             // calculate flux through TRIANGLE
                        }
                }

                for(i=0;i<N_POINTS;i++){
                        for(j=0;j<N_POINTS;j++){                                     // loop over all vertices
                                POINTS[i][j].update_u_variables();                   // update the fluid state at vertex
                                POINTS[i][j].con_to_prim();                          // convert these to their corresponding conserved
                                POINTS[i][j].recalculate_pressure();                 // caclulate pressure from new conserved values
                                //POINTS[i][j].prim_to_con();
                                POINTS[i][j].reset_du();                             // reset du value to zero for next timestep
                                //POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);        // calculate next timestep based on new state
                                //if(POSSIBLE_DT<NEXT_DT){NEXT_DT = POSSIBLE_DT;}
                        }
                }

                T+=DT;                                                            // increment time
                l+=1;                                                             // increment step number

        }

        output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DX);      // write out final state
        close_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);

        return 0;
}
