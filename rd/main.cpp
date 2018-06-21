#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>

#include "cblas.h"
#include "lapacke.h"
#include "inverse.cpp"

#include "constants.h"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#include "setup2D.cpp"
#include "io2D.cpp"
#endif

int main(){

        int i, j, l = 0;                                          // ******* decalare varaibles and std::vectors ******
        double DX, DY, DT, T = 0.0;                                // DX = space step,DT = timestep,t = time,CFL = CFL condition
        double NEXT_TIME  = 0.0;                                   // T_TOT = total time,NEXT_TIME = time of next snapshot
        VERTEX                               NEW_VERTEX;
        TRIANGLE                             NEW_TRIANGLE;
        std::vector<VERTEX>                  X_POINTS;
        std::vector<std::vector<VERTEX> >    POINTS;
        std::vector<TRIANGLE>                X_MESH;
        std::vector<std::vector<TRIANGLE> >  MESH;
        std::vector<VERTEX>::iterator        IT_VERT;
        std::vector<TRIANGLE>::iterator      IT_TRIANGLE;
        double NEXT_DT = T_TOT, POSSIBLE_DT;               // TOTAL_DENSITY = total density in box

        /****** Setup initial conditions of one dimensional tube ******/

#ifdef LDA_SCHEME
        std::cout << "Using LDA Scheme" << std::endl;
#endif

#ifdef N_SCHEME
        std::cout << "Using N Scheme" << std::endl;
#endif

        std::cout << "Building grid of vertices" << std::endl;

#ifdef TWO_D
        for(j=0; j<N_POINTS_Y; j++){
                for(i=0; i<N_POINTS_X; i++){
                        NEW_VERTEX = setup_vertex(i,j,DX,DY);               // call VERTEX setup routine
                        X_POINTS.push_back(NEW_VERTEX);                              // add new VERTEX to std::vector of vertices in this row
                        X_POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);                // check dt is min required by CFL
                        if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                }
                POINTS.push_back(X_POINTS);
                X_POINTS.clear();
        }
#endif

        /****** Setup MESH ******/

         std::cout << "Assigning vertices to triangles" << std::endl;

        for(j=0; j<2*N_POINTS_Y; j++){
                for(i=0; i<N_POINTS_X; i++){
                        NEW_TRIANGLE = setup_triangle(i,j,POINTS);
                        X_MESH.push_back(NEW_TRIANGLE);
                }
                MESH.push_back(X_MESH);
                X_MESH.clear();
        }



        std::cout << "Finished triangle setup" << std::endl;

        /****** Loop over time until total time T_TOT is reached ******/

        std::ofstream DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE;

        open_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);               // open output files

        std::cout << "Mesh Size =\t" << MESH[0].size() << '\t' << MESH.size() << std::endl;
        std::cout << "Evolving fluid ..." << std::endl;

        while(T<T_TOT){

                DT = NEXT_DT;                                                     // set timestep based oncaclulation from previous timestep

#ifdef FIXED_DT
                DT = 0.00001;
#endif

                std::cout << "STEP =\t" << l << "\tTIME =\t" << T << "\tTIMESTEP =\t" << DT << std::endl;

                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DY);
                }

#ifdef DEBUG
                // std::cout << std::fixed;
                // std::cout << std::setprecision(10);
                std::cout << "Calculating first half time step change" << std::endl;
#endif

                
                for(j=0;j<2*N_POINTS_Y;j++){                                        // loop over all triangles in MESH
                        for(i=0;i<N_POINTS_X;i++){ 
                                MESH[j][i].calculate_first_half(T, DT, DX, DY);             // calculate flux through TRIANGLE
                        }
                }

                for(j=0;j<N_POINTS_Y;j++){
                        for(i=0;i<N_POINTS_X;i++){                                     // loop over all vertices
                                POINTS[j][i].update_u_half();                        // update the half time state
                                POINTS[j][i].con_to_prim_half();
                                POINTS[j][i].recalculate_pressure_half();            // caclulate pressure from new conserved values
                                POINTS[j][i].reset_du_half();                             // reset du value to zero for next timestep
                        }
                }

                for(j=0;j<2*(N_POINTS_Y);j++){                                         // loop over all triangles in MESH
                        for(i=0;i<N_POINTS_X;i++){
                                MESH[j][i].calculate_second_half(T, DT, DX, DY);             // calculate flux through TRIANGLE
                        }
                }

                NEXT_DT = T_TOT - (T + DT);        // set next timestep to max possible value (time remaining to end)

                for(j=0;j<N_POINTS_Y;j++){
                        for(i=0;i<N_POINTS_X;i++){                                     // loop over all vertices
                                POINTS[j][i].update_u_variables();                   // update the fluid state at vertex
                                POINTS[j][i].con_to_prim();                          // convert these to their corresponding conserved
                                POINTS[j][i].recalculate_pressure();                 // caclulate pressure from new conserved values
                                POINTS[j][i].reset_du();                             // reset du value to zero for next timestep
                                POINTS[j][i].calc_next_dt(DX,CFL,POSSIBLE_DT);        // calculate next timestep based on new state
                                if(POSSIBLE_DT<NEXT_DT){NEXT_DT = POSSIBLE_DT;}
                        }
                }

                T+=DT;                                                            // increment time
                l+=1;                                                             // increment step number

        }

        output_state(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE, POINTS, T, DT, DX, DX);      // write out final state
        close_files(DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, DU_FILE);

        return 0;
}
