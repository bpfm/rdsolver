#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <omp.h> 

#include "constants.h"

#include "cblas.h"
#include "lapacke.h"
#include "inverse.cpp"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#include "setup2D.cpp"
#include "io2D.cpp"
#include "gravity2D.cpp"
#endif

int main(){
        /*
        Setup and run simulation from input file constants.h, using precalculated triangulation
        */

        int i, j, l = 0, m;                                           // ******* decalare varaibles and vectors ******
        double DX, DY, DT, T = 0.0;                                // DX           = space step,DT = timestep,t = time
        double NEXT_TIME = 0.0;                                    // NEXT_TIME    = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT = T_TOT;                       // NEXT_DT.     = timestep for upcoming time iteration
        double MIN_DT;
        VERTEX                               NEW_VERTEX;           // NEW_VERTEX   = dummy variable for setting up vertices
        TRIANGLE                             NEW_TRIANGLE;         // NEW_TRIABLE  = dummy variable for setting up triangles
        std::vector<VERTEX>                  RAND_POINTS;             // X_POINTS     = vector of x vertices
        std::vector<TRIANGLE>                RAND_MESH;               // X_MESH       = vector of x triangles
        double SNAP_ID = 0;


        // Initialise seed for random number generator (rand)
        std::srand(68315);

        /****** Setup initial conditions of one dimensional tube ******/

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
                NEW_VERTEX.reset_len_vel_sum();
                RAND_POINTS.push_back(NEW_VERTEX);
        }


        /****** Setup mesh ******/

        std::cout << "Reading QHULL triangles ..." << std::endl;

        std::cout << "Number of triangles = " << N_TRIANG << std::endl;

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = qhull_read_triangles_line(TRIANGLES_FILE,RAND_POINTS);
                NEW_TRIANGLE.set_tbin(1);
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

        CGAL_FILE_NAME = "Delaunay2D.txt";

        CGAL_FILE.open(CGAL_FILE_NAME);

        N_POINTS = cgal_read_positions_header(CGAL_FILE);

        std::cout << "Number of vertices = " << N_POINTS << std::endl;

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = cgal_read_positions_line(CGAL_FILE);
                NEW_VERTEX.reset_len_vel_sum();
                NEW_VERTEX.set_id(i);
                RAND_POINTS.push_back(NEW_VERTEX);
        }

        /****** Setup mesh ******/

        std::cout << "Reading CGAL triangles ..." << std::endl;

        N_TRIANG = cgal_read_triangles_header(CGAL_FILE);

        std::cout << "Number of triangles = " << N_TRIANG << std::endl;

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = cgal_read_triangles_line(CGAL_FILE,RAND_POINTS,j);
                NEW_TRIANGLE.set_tbin(1);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

#endif
#endif

#ifdef SEDOV
        double ETOT = 0.0,ETOT_AIM = 300000.0,PRESSURE_AIM;
        for(i=0; i<N_POINTS; ++i){
                if((RAND_POINTS[i].get_x()-5.0)*(RAND_POINTS[i].get_x()-5.0) + (RAND_POINTS[i].get_y()-5.0)*(RAND_POINTS[i].get_y()-5.0) < R_BLAST*R_BLAST){
                        AREA_CHECK = AREA_CHECK + RAND_POINTS[i].get_dual();
                }
        }
        for(i=0; i<N_POINTS; ++i){
                if((RAND_POINTS[i].get_x()-5.0)*(RAND_POINTS[i].get_x()-5.0) + (RAND_POINTS[i].get_y()-5.0)*(RAND_POINTS[i].get_y()-5.0) < R_BLAST*R_BLAST){
                        PRESSURE_AIM = (ETOT_AIM * GAMMA_1 / RAND_POINTS[i].get_dual()) * (RAND_POINTS[i].get_dual() / (AREA_CHECK));
                        RAND_POINTS[i].set_pressure(PRESSURE_AIM);
                        ETOT = ETOT + RAND_POINTS[i].get_pressure()*RAND_POINTS[i].get_dual()/GAMMA_1;
                        std::cout << POINT_CHECK << "\t" << PRESSURE_AIM << "\t" << RAND_POINTS[i].get_pressure() << "\t" << ETOT << std::endl;
                        RAND_POINTS[i].setup_specific_energy();
                        RAND_POINTS[i].prim_to_con();
                }
        }
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

        std::cout << std::fixed;
        std::cout << std::setprecision(6);

        int TBIN, TBIN_CURRENT = 0;

        NEXT_DT = 0.0;                                                            // set first timestep to zero

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

#ifdef PARA_RES
                #pragma omp parallel for
#endif
                for(j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        TBIN = RAND_MESH[j].get_tbin();
                        if(TBIN_CURRENT == 0 or (TBIN_CURRENT == 1 and  TBIN == 1) \
                                             or (TBIN_CURRENT == 2 and (TBIN == 2 or TBIN == 1)) \
                                             or (TBIN_CURRENT == 3 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 4 and (TBIN == 4 or TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 5 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 6 and (TBIN == 2 or TBIN == 1)) \
                                             or (TBIN_CURRENT == 7 and  TBIN == 1)\
                                             ){
                                // std::cout << TBIN_CURRENT << "\t" << RAND_MESH[j].get_tbin() <<std::endl;
                                RAND_MESH[j].calculate_first_half(T);
                        }
                        // RAND_MESH[j].calculate_first_half(T);                                                 // calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update_half(DT);
                }


#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_half();                        // update the half time state
                        RAND_POINTS[i].reset_du_half();                        // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values_half();
                        RAND_POINTS[i].con_to_prim_half();
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
                        RAND_POINTS[i].reset_du();                             // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values();
                        RAND_POINTS[i].con_to_prim();                          // convert these to their corresponding conserved
                }

#ifdef SELF_GRAVITY
                direct_gravity(RAND_POINTS, N_POINTS, DT);
#endif

                if(TBIN_CURRENT == 0){
                        for(j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                                RAND_MESH[j].calculate_len_vel_contribution();         // calculate flux through TRIANGLE
                        }
                        NEXT_DT = T_TOT - (T + DT);        // set next timestep to max possible value (time remaining to end)å
                        for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                                POSSIBLE_DT = RAND_POINTS[i].calc_next_dt();           // calculate next timestep based on new state
                                if(POSSIBLE_DT < NEXT_DT){NEXT_DT = POSSIBLE_DT;}
                                RAND_POINTS[i].reset_len_vel_sum();
                                // RAND_POINTS[i].set_tbin_local(N_TBINS);
                        }
                        for(j=0;j<N_TRIANG;++j){                                        // bin triangles by minimum timestep of vertices
                                MIN_DT = RAND_MESH[j].get_vertex_0()->get_dt_req();
                                if(RAND_MESH[j].get_vertex_1()->get_dt_req() < MIN_DT){MIN_DT = RAND_MESH[j].get_vertex_1()->get_dt_req();}
                                if(RAND_MESH[j].get_vertex_2()->get_dt_req() < MIN_DT){MIN_DT = RAND_MESH[j].get_vertex_2()->get_dt_req();}
                                if(MIN_DT < 2.0*NEXT_DT){
                                        RAND_MESH[j].set_tbin(1);
                                        // std::cout << 1 << std::endl;
                                }else if(MIN_DT > 2.0*NEXT_DT and MIN_DT < 4.0*NEXT_DT){
                                        RAND_MESH[j].set_tbin(2);
                                        // std::cout << 2 << std::endl;
                                }else if(MIN_DT > 4.0*NEXT_DT and MIN_DT < 8.0*NEXT_DT){
                                        RAND_MESH[j].set_tbin(4);
                                        // std::cout << 4 << std::endl;
                                }else{
                                        RAND_MESH[j].set_tbin(8);
                                        // std::cout << 8 << std::endl
                                }
                                // RAND_MESH[j].send_tbin_limit();
                        }
                        // for(j=0;j<N_TRIANG;++j){
                        //         RAND_MESH[j].check_tbin();
                        // }
                }
                // std::cout << NEXT_DT << std::endl;
                // std::cout << TBIN_CURRENT << std::endl;

                TBIN_CURRENT = (TBIN_CURRENT + 1) % N_TBINS;                     // increment time step bin
                T += DT;                                                         // increment time
                l += 1;                                                          // increment step number
        }

        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);

        return 0;
}
