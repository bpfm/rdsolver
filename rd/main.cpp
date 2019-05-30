#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
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

        int i, j, l = 0;                                           // ******* decalare varaibles and vectors ******
        double DX, DY, DT, T = 0.0;                                // DX           = space step,DT = timestep,t = time
        double NEXT_TIME = 0.0;                                    // NEXT_TIME    = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT = T_TOT;                       // NEXT_DT.     = timestep for upcoming time iteration
        VERTEX                               NEW_VERTEX;           // NEW_VERTEX   = dummy variable for setting up vertices
        TRIANGLE                             NEW_TRIANGLE;         // NEW_TRIABLE  = dummy variable for setting up triangles
        std::vector<VERTEX>                  RAND_POINTS;             // X_POINTS     = vector of x vertices
        std::vector<TRIANGLE>                RAND_MESH;               // X_MESH       = vector of x triangles


        // Initialise seed for random number generator (rand)
        std::srand(68315);

        /****** Setup initial conditions of one dimensional tube ******/

        std::cout << std::fixed;
        std::cout << std::setprecision(9);

        std::cout << "*********************************************************" << std::endl;

#ifdef LDA_SCHEME
        std::cout << "Using LDA Scheme" << std::endl;
#endif

#ifdef N_SCHEME
        std::cout << "Using N Scheme" << std::endl;
#endif

#ifdef BLENDED
        std::cout << "Using B Scheme" << std::endl;
#endif

#ifdef FIRST_ORDER
        std::cout << "Using 1st order" << std::endl;
#else
        std::cout << "Using 2nd order" << std::endl;
#endif

        std::cout << "Building grid of vertices" << std::endl;

        /****** Setup Vertices ******/

#ifdef READ_IC
#ifdef QHULL_IC
        int N_POINTS, N_TRIANG;
        std::string   POSITIONS_FILE_NAME, TRIANGLES_FILE_NAME;
        std::ifstream POSITIONS_FILE, TRIANGLES_FILE;

        /****** Setup vertices ******/

        std::cout << "Reading QHULL vertex positions ..." << std::endl;

        POSITIONS_FILE_NAME = "triangulation/points.txt";
        TRIANGLES_FILE_NAME = "triangulation/ordered_triangles.txt";

        POSITIONS_FILE.open(POSITIONS_FILE_NAME);
        TRIANGLES_FILE.open(TRIANGLES_FILE_NAME);

        N_POINTS = qhull_read_positions_header(POSITIONS_FILE);
        N_TRIANG = qhull_read_triangles_header(TRIANGLES_FILE);

        std::cout << "Number of vertices = " << N_POINTS << std::endl;

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = qhull_read_positions_line(POSITIONS_FILE);
                RAND_POINTS.push_back(NEW_VERTEX);
        }


        /****** Setup mesh ******/

        std::cout << "Reading QHULL triangles ..." << std::endl;

        std::cout << "Number of triangles = " << N_TRIANG << std::endl;

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = qhull_read_triangles_line(TRIANGLES_FILE,RAND_POINTS);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

        POSITIONS_FILE.close();
        TRIANGLES_FILE.close();
#endif
#ifdef CGAL_IC
        int N_POINTS, N_TRIANG;
        std::string   CGAL_FILE_NAME;
        std::ifstream CGAL_FILE;

        /****** Setup vertices ******/

        std::cout << "Reading CGAL vertex positions ..." << std::endl;

        CGAL_FILE_NAME = "triangulation/cgal/output.txt";

        CGAL_FILE.open(CGAL_FILE_NAME);

        N_POINTS = cgal_read_positions_header(CGAL_FILE);

        std::cout << "Number of vertices = " << N_POINTS << std::endl;

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = cgal_read_positions_line(CGAL_FILE);
                RAND_POINTS.push_back(NEW_VERTEX);
        }

        /****** Setup mesh ******/

        std::cout << "Reading CGAL triangles ..." << std::endl;

        N_TRIANG = cgal_read_triangles_header(CGAL_FILE);

        std::cout << "Number of triangles = " << N_TRIANG << std::endl;

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = cgal_read_triangles_line(CGAL_FILE,RAND_POINTS);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

#endif        
#endif

        /****** Set initial timestep  ******/

        std::cout << "Finding initial timestep ..." << std::endl;

        for(i=0; i<N_POINTS; ++i){
                NEXT_DT = RAND_POINTS[i].calc_next_dt();      // check dt is min required by CFL
                if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
        }

        std::ofstream POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, TEMP;

        // TEMP.open("output/temp.txt");

        open_files(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN);               // open output files

        std::cout << "Checking mesh size ..." << std::endl;

        std::cout << "Mesh Size =\t" << RAND_MESH.size() << std::endl;

        std::cout << "Evolving fluid ..." << std::endl;

        /****** Loop over time until total time T_TOT is reached ******/

        while(T<T_TOT){

                DT = NEXT_DT;                                                     // set timestep based oncaclulation from previous timestep

#ifdef FIXED_DT
                DT = 0.00001;
#endif

                std::cout << "STEP =\t" << l << "\tTIME =\t" << T << "\tTIMESTEP =\t" << DT << std::endl;

                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        output_state(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, RAND_POINTS, T, DT, N_POINTS);
                }

// #ifdef DEBUG
                // std::cout << std::fixed;
                // std::cout << std::setprecision(6);
                // std::cout << "Calculating first half time step change" << std::endl;
// #endif

                for(j=0;j<N_TRIANG;++j){                                        // loop over all triangles in MESH
                        RAND_MESH[j].calculate_first_half(T, DT);             // calculate flux through TRIANGLE

                }

                for(i=0;i<N_POINTS;++i){                                   // loop over all vertices
                        RAND_POINTS[i].update_u_half();                        // update the half time state
                        RAND_POINTS[i].con_to_prim_half();
                        RAND_POINTS[i].reset_du_half();                             // reset du value to zero for next timestep

                }

                for(j=0;j<N_TRIANG;++j){                                         // loop over all triangles in MESH
                        RAND_MESH[j].calculate_second_half(T, DT);             // calculate flux through TRIANGLE
                }

                NEXT_DT = T_TOT - (T + DT);        // set next timestep to max possible value (time remaining to end)

                for(i=0;i<N_POINTS;++i){                                   // loop over all vertices
                        RAND_POINTS[i].update_u_variables();                   // update the fluid state at vertex
                        RAND_POINTS[i].con_to_prim();                          // convert these to their corresponding conserved
                        RAND_POINTS[i].reset_du();                             // reset du value to zero for next timestep
                        POSSIBLE_DT = RAND_POINTS[i].calc_next_dt();        // calculate next timestep based on new state
                        if(POSSIBLE_DT<NEXT_DT){NEXT_DT = POSSIBLE_DT;}
                }

                T+=DT;                                                            // increment time
                l+=1;                                                             // increment step number

        }

        output_state(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN, RAND_POINTS, T, DT, N_POINTS);      // write out final state
        close_files(POSITIONS, DENSITY_MAP, PRESSURE_MAP, VELOCITY_MAP, CENTRAL_COLUMN);

        return 0;
}
