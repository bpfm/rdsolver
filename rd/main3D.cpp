#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <omp.h> 


#include "constants3D.h"

#include "cblas.h"
#include "lapacke.h"
#include "inverse.cpp"

#ifdef THREE_D
#include "vertex3D.h"
#include "triangle3D.h"
#include "setup3D.cpp"
#include "io3D.cpp"
// #include "gravity3D.cpp"
#endif

int main(){
        /*
        Setup and run simulation from input file constants.h, using precalculated triangulation
        */

        int i, j, k, l = 0, m;                                     // ******* decalare varaibles and vectors ******
        double DT, T = 0.0;                                        // DT = timestep,t = time
        double NEXT_TIME = 0.0;                                    // NEXT_TIME    = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT = T_TOT;               // NEXT_DT.     = timestep for upcoming time iteration
        double MIN_DT;
        VERTEX                               NEW_VERTEX;           // NEW_VERTEX   = dummy variable for setting up vertices
        TRIANGLE                             NEW_TRIANGLE;         // NEW_TRIABLE  = dummy variable for setting up triangles
        std::vector<VERTEX>                  RAND_POINTS;          // X_POINTS     = vector of x vertices
        std::vector<TRIANGLE>                RAND_MESH;            // X_MESH       = vector of x triangles
        double SNAP_ID = 0;


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

        std::ofstream LOGFILE;
        LOGFILE << std::setprecision(12);
        LOGFILE.open("output/log.txt");

        /****** Setup Vertices ******/

#ifdef READ_IC
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
                NEW_VERTEX.reset_len_vel_sum();
                RAND_POINTS.push_back(NEW_VERTEX);
        }

        /****** Setup mesh ******/

        std::cout << "Reading CGAL triangles ..." << std::endl;

        N_TRIANG = cgal_read_triangles_header(CGAL_FILE);

        std::cout << "Number of triangles = " << N_TRIANG << std::endl;

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = cgal_read_triangles_line(CGAL_FILE,RAND_POINTS);
                NEW_TRIANGLE.set_tbin(1);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

#endif
#endif

        /****** Set initial timestep  ******/

        std::cout << "Finding initial timestep ..." << std::endl;

        for(j=0;j<N_TRIANG;++j){                                           // loop over all triangles in MESH
                RAND_MESH[j].calculate_len_vel_contribution();             // calculate flux through TRIANGLE
        }

        for(i=0; i<N_POINTS; ++i){
                NEXT_DT = RAND_POINTS[i].calc_next_dt();      // check dt is min required by CFL
                if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                RAND_POINTS[i].reset_len_vel_sum();
        }

        std::cout << "Checking mesh size ..." << std::endl;
        std::cout << "Mesh Size =\t" << RAND_MESH.size() << std::endl;
        std::cout << "Evolving fluid ..." << std::endl;

        /****** Loop over time until total time T_TOT is reached *****************************************************************************************************/

        int TBIN_CURRENT = 0;
        int TBIN;
        int ACTIVE, ACTIVE_ID = 0;

        NEXT_DT = 0.0;

        while(T<T_TOT){

                DT = NEXT_DT;                                                     // set timestep based oncaclulation from previous timestep

#ifdef FIXED_DT
                DT = DT_FIX;
#endif

                std::cout << "STEP =\t" << l << "\tTIME =\t" << T << "\tTIMESTEP =\t" << DT << "\t" << 100.0*T/T_TOT << " %" <<  "\r" << std::flush;

                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);
                        write_active(RAND_MESH, N_TRIANG, SNAP_ID, TBIN_CURRENT);
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        SNAP_ID ++;
                }

// #ifdef DEBUG
                // std::cout << std::fixed;
                // std::cout << std::setprecision(6);
                // std::cout << "Calculating first half time step change" << std::endl;
// #endif
                ACTIVE = 0;

#ifdef PARA_RES
                #pragma omp parallel for
#endif
                for(j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        RAND_MESH[j].calculate_first_half(T);                                                 // calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update_half(DT);
                }


#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_half();                        // update the half time state
                        RAND_POINTS[i].con_to_prim_half();
                        RAND_POINTS[i].reset_du_half();                        // reset du value to zero for next timestep
                }

#ifdef PARA_RES
                #pragma omp parallel for
#endif
                for(j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                        RAND_MESH[j].calculate_second_half(T, DT);             // calculate flux through TRIANGLE
                }

#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_variables();                   // update the fluid state at vertex
                        RAND_POINTS[i].con_to_prim();                          // convert these to their corresponding conserved
                        RAND_POINTS[i].reset_du();                             // reset du value to zero for next timestep
                }

// #ifdef SELF_GRAVITY
//                 direct_gravity(RAND_POINTS, N_POINTS, DT);
// #endif

                for(j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                        RAND_MESH[j].calculate_len_vel_contribution();         // calculate flux through TRIANGLE
                }

                T += DT;                                                         // increment time
                l += 1;                                                          // increment step number
        }

        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);

        return 0;
}
