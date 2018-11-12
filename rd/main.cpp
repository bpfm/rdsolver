#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>

#include "constants.h"

#include "cblas.h"
#include "lapacke.h"
#include "inverse.cpp"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#include "setup2D.cpp"
#include "io2D.cpp"
#endif

int main(){

        int i, j, l = 0;                                          // ******* decalare varaibles and std::vectors ******
        double DX, DY, DT, T = 0.0;                                // DX = space step,DT = timestep,t = time
        double NEXT_TIME = 0.0;                                   // NEXT_TIME = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT;                       // NEXT_DT = timestep for upcoming time iteration
        VERTEX                               NEW_VERTEX;
        TRIANGLE                             NEW_TRIANGLE;
        std::vector<VERTEX>                  X_POINTS;
        std::vector<std::vector<VERTEX> >    POINTS;
        std::vector<TRIANGLE>                X_MESH;
        std::vector<std::vector<TRIANGLE> >  MESH;
        std::vector<VERTEX>::iterator        IT_VERT;
        std::vector<TRIANGLE>::iterator      IT_TRIANGLE;


        /****** Setup initial conditions of one dimensional tube ******/

#ifdef LDA_SCHEME
        std::cout << "Using LDA Scheme" << std::endl;
#endif

#ifdef N_SCHEME
        std::cout << "Using N Scheme" << std::endl;
#endif

        std::cout << "Building grid of vertices" << std::endl;

#ifdef GENERATE_IC
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
#endif

#ifdef READ_IC
        std::string   IC_FILE_NAME;
        std::ifstream IC_FILE;

        IC_FILE_NAME = "ic.txt";

        IC_FILE.open(IC_FILE_NAME);

        for(j=0; j<N_POINTS_Y; j++){
                for(i=0; i<N_POINTS_X; i++){
                        NEW_VERTEX = READ_IC_LINE(IC_FILE);
                        X_POINTS.push_back(NEW_VERTEX);                     // add new VERTEX to std::vector of vertices in this row
                        X_POINTS[i].calc_next_dt(DX,CFL,POSSIBLE_DT);       // check dt is min required by CFL
                        if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                }
                POINTS.push_back(X_POINTS);
                X_POINTS.clear();
        }

        IC_FILE.close();
#endif

        /****** Setup MESH ******/

#ifdef MESH_TEST
        std::ofstream MESH_POS;
        MESH_POS.open("mesh_pos.txt");
#endif

         std::cout << "Assigning vertices to triangles" << std::endl;

        for(j=0; j<2*N_POINTS_Y; j++){
                for(i=0; i<N_POINTS_X; i++){
                        NEW_TRIANGLE = setup_triangle(i,j,POINTS);
                        X_MESH.push_back(NEW_TRIANGLE);
#ifdef MESH_TEST
                        // if(i>0 and i<N_POINTS_X and j>0 and j<2*N_POINTS_Y){
                                // std::cout << i << "\t" << j << std::endl;
                                double X_CENTRE = (NEW_TRIANGLE.get_vertex_0()->get_x() + NEW_TRIANGLE.get_vertex_1()->get_x() + NEW_TRIANGLE.get_vertex_2()->get_x())/3.0;
                                double Y_CENTRE = (NEW_TRIANGLE.get_vertex_0()->get_y() + NEW_TRIANGLE.get_vertex_1()->get_y() + NEW_TRIANGLE.get_vertex_2()->get_y())/3.0;
                                double X_DIST_0 = X_CENTRE - NEW_TRIANGLE.get_vertex_0()->get_x();
                                double X_DIST_1 = X_CENTRE - NEW_TRIANGLE.get_vertex_1()->get_x();
                                double X_DIST_2 = X_CENTRE - NEW_TRIANGLE.get_vertex_2()->get_x();
                                double Y_DIST_0 = Y_CENTRE - NEW_TRIANGLE.get_vertex_0()->get_y();
                                double Y_DIST_1 = Y_CENTRE - NEW_TRIANGLE.get_vertex_1()->get_y();
                                double Y_DIST_2 = Y_CENTRE - NEW_TRIANGLE.get_vertex_2()->get_y();
                                double LEN_0 = sqrt(X_DIST_0*X_DIST_0 + Y_DIST_0*Y_DIST_0);
                                double LEN_1 = sqrt(X_DIST_1*X_DIST_1 + Y_DIST_0*Y_DIST_0);
                                double LEN_2 = sqrt(X_DIST_2*X_DIST_2 + Y_DIST_0*Y_DIST_0);
                                MESH_POS << NEW_TRIANGLE.get_vertex_0()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_0()->get_y() << "\t" << NEW_TRIANGLE.get_vertex_1()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_1()->get_y() << "\t" << NEW_TRIANGLE.get_vertex_2()->get_x() << "\t" << NEW_TRIANGLE.get_vertex_2()->get_y() << "\t" << X_CENTRE << "\t" << Y_CENTRE << "\t" << LEN_0 << "\t" << LEN_1 << "\t" << LEN_2 << std::endl;
                        // }
#endif
                }
                MESH.push_back(X_MESH);
                X_MESH.clear();
        }

#ifdef MESH_TEST
        exit(0);
#endif

        std::cout << "Finished triangle setup" << std::endl;

        /****** Loop over time until total time T_TOT is reached ******/

        std::ofstream POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, GENERATED_IC;

        open_files(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, GENERATED_IC);               // open output files

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
                        output_state(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, GENERATED_IC, POINTS, T, DT, DX, DY);
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
#ifdef UPDATE_TEST
                                POINTS[j][i].output_update_count();
#endif
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

#ifdef UPDATE_TEST
                exit(0);
#endif
        }

        output_state(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, GENERATED_IC, POINTS, T, DT, DX, DX);      // write out final state
        close_files(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, GENERATED_IC);

        return 0;
}
