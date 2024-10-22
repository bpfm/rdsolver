/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <exception>
#include "constants.h"
#include "all_functions.h"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#endif

#ifdef SEDOV
    int POINT_CHECK = 0;
#endif

double GLOBAL_PRESSURE_MAX, GLOBAL_PRESSURE_MIN = 0.0;
double GLOBAL_VELOCITY_SCALE;
double MESH_RESOLUTION = 0.0;

//using namespace std;
namespace fs = std::filesystem;

int main(int ARGC, char *ARGV[]){
        /*
        Setup and run simulation from input file constants.h, using precalculated triangulation
        */
        auto program_start = std::chrono::system_clock::now();
        std::time_t program_start_time = std::chrono::system_clock::to_time_t(program_start);
        std::cout<<"start computation at "<<std::ctime(&program_start_time);

        int i, j, l = 0, m;                                        // ******* declare variables and vectors ******
        int SNAP_ID = 0;
        double DT, T = 0.0;                                        //
        double NEXT_TIME = 0.0;                                    // NEXT_TIME    = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT = T_TOT;               // NEXT_DT.     = timestep for upcoming time iteration
        double MIN_DT;
        VERTEX                               NEW_VERTEX;           // NEW_VERTEX   = dummy variable for setting up vertices
        TRIANGLE                             NEW_TRIANGLE;         // NEW_TRIABLE  = dummy variable for setting up triangles
        std::vector<VERTEX>                  RAND_POINTS;          // X_POINTS     = vector of x vertices
        std::vector<TRIANGLE>                RAND_MESH;            // X_MESH       = vector of x triangles

        // Initialise seed for random number generator (rand)
        std::srand(68315);

        // read_parameter_file(ARGC, ARGV);

        /****** Setup simulation options ******/

        printf("*********************************************************\n");

        printf("LAIRDS 2D\n");

#ifdef LDA_SCHEME
        printf("Using LDA Scheme\n");
#endif

#ifdef N_SCHEME
        printf("Using N Scheme\n");
#endif

#ifdef BLENDED
        printf("Using B Scheme\n");
#endif

#ifdef BLENDED_X
    printf("Using Bx Scheme\n");
#endif


#ifdef FIRST_ORDER
        printf("Using 1st order\n");
#else
        printf("Using 2nd order\n");
#endif

        printf("Building vertices and mesh\n");


#ifdef PARA_RES
        omp_set_dynamic(0);
        omp_set_num_threads(8);
#pragma omp parallel
        {
            int thread_ID = omp_get_num_threads();
            cout<<"number of threads: "<<thread_ID<<endl;
        }
#endif


        fs::path sourceFile = "constants.h";
        fs::path targetParent = OUT_DIR;
        auto target = targetParent / sourceFile.filename();
        try // If you want to avoid exception handling, then use the error code overload of the following functions.
        {
            fs::create_directories(targetParent); // Recursively create target directory if not existing.
            fs::copy_file(sourceFile, target, fs::copy_options::overwrite_existing);
        }
        catch (std::exception& e) // Not using fs::filesystem_error since std::bad_alloc can throw too.
        {
            std::cout << e.what();
        }

        std::ofstream LOGFILE;
        LOGFILE << std::setprecision(12);
        LOGFILE.open(LOG_DIR);

#ifdef READ_IC
#ifdef QHULL_IC
        int N_POINTS, N_TRIANG;
        std::string   POSITIONS_FILE_NAME, TRIANGLES_FILE_NAME;
        std::ifstream POSITIONS_FILE, TRIANGLES_FILE;

        /****** Setup vertices ******/

        printf("Reading QHULL vertex positions ...\n");

        POSITIONS_FILE_NAME = "triangulation/points.txt";
        TRIANGLES_FILE_NAME = "triangulation/ordered_triangles.txt";

        POSITIONS_FILE.open(POSITIONS_FILE_NAME);
        TRIANGLES_FILE.open(TRIANGLES_FILE_NAME);

        N_POINTS = qhull_read_positions_header(POSITIONS_FILE);
        N_TRIANG = qhull_read_triangles_header(TRIANGLES_FILE);

        printf("Number of vertices = %d\n", N_POINTS);

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = qhull_read_positions_line(POSITIONS_FILE);
                NEW_VERTEX.reset_len_vel_sum();
                RAND_POINTS.push_back(NEW_VERTEX);
        }


        /****** Setup mesh ******/

        printf("Reading QHULL triangles ...");

        printf("Number of triangles = %d\n", N_TRIANG);

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
        std::ifstream CGAL_FILE;

        /****** Setup vertices ******/


        printf("Reading CGAL vertex positions ...");

        CGAL_FILE.open(CGAL_FILE_NAME);
        N_POINTS = cgal_read_positions_header(CGAL_FILE);

        printf("Number of vertices = %d\n", N_POINTS);

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = cgal_read_positions_line(CGAL_FILE);
                NEW_VERTEX.reset_len_vel_sum();
                NEW_VERTEX.set_id(i);
                RAND_POINTS.push_back(NEW_VERTEX);
        }

        /****** Setup mesh ******/

        printf("Reading CGAL triangles ...");

        N_TRIANG = cgal_read_triangles_header(CGAL_FILE);

        printf("Number of triangles = %d\n", N_TRIANG);

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = cgal_read_triangles_line(CGAL_FILE,RAND_POINTS,j);
                NEW_TRIANGLE.set_tbin(1);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

#endif
#endif

#ifdef SEDOV
        /****** Inject pressure for Sedov test  ******/

        double AREA_CHECK  = 0.0;
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
                        printf("%f\t%f\t%f\n", PRESSURE_AIM, RAND_POINTS[i].get_pressure(), ETOT);
                        RAND_POINTS[i].setup_specific_energy();
                        RAND_POINTS[i].prim_to_con();
                }
        }
#endif

        /****** Set initial timestep  ******/

        printf("Finding initial timestep ...");

        for(j=0;j<N_TRIANG;++j){                                           // loop over all triangles in MESH
                RAND_MESH[j].calculate_len_vel_contribution();             // calculate flux through TRIANGLE
        }

        for(i=0; i<N_POINTS; ++i){
                NEXT_DT = RAND_POINTS[i].calc_next_dt();      // check dt is min required by CFL
                if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                RAND_POINTS[i].reset_len_vel_sum();
        }

        printf("Checking mesh size ...");
        printf("Mesh Size = %d\n",int(RAND_MESH.size()));
        printf("Evolving fluid ...\n");

        int TBIN, TBIN_CURRENT = 0;
        //NEXT_DT = 0.0;                                                            // set first timestep to zero

        /****** Loop over time until total time T_TOT is reached *****************************************************************************************************/
        while(T<T_TOT){

        /****** Update time step to new value ******/
                DT = NEXT_DT;                                                     // set timestep based on caclulation from previous timestep

#ifdef FIXED_DT
        /****** Reset time step if fixed ******/
                DT = DT_FIX;
#endif
                printf("STEP =\t%d\tTIME =\t%f\tTIMESTEP =\t%f\t%f/100\n", l, T, DT, 100.0*T/T_TOT);

            /****** Write snapshot *****************************************************************************************************/
                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);
                        //write_active(RAND_MESH, N_TRIANG, SNAP_ID, TBIN_CURRENT);
#if defined(BLENDED) or defined(BLENDED_X)
                        write_vertex_list(RAND_POINTS, T, N_POINTS,SNAP_ID);
                        write_blending_coeff(RAND_MESH, T, N_TRIANG, SNAP_ID);
#endif
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        SNAP_ID ++;
                }

        /****** 1st order update ***************************************************************************************************/
#ifdef BLENDED_X
                set_shock_sensor(N_POINTS,RAND_POINTS);
#endif


#ifdef DRIFT
                /****** Update residual for active bins (Drift method) ******/
                drift_update_half(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH);
#endif
#ifdef JUMP
                /****** Update residual for active bins (Jump method) ******/
                jump_update_half(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH);
#endif




#if !defined(DRIFT) && !defined(JUMP)
#ifdef PARA_RES
                #pragma omp parallel for
#endif
                /****** Update residual for all bins (No adaptive method) ******/
                for(j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        RAND_MESH[j].calculate_first_half(T,DT);                                                 // calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update_half();
                }
#endif

#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_half();                        // update the half time state
                        RAND_POINTS[i].reset_du_half();                        // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values_half();
                        RAND_POINTS[i].con_to_prim_half();
                }

        /****** 2nd order update ***************************************************************************************************/

#ifdef DRIFT
                /****** Update residual for active bins (Drift method) ******/
                drift_update(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH);
#endif

#if !defined(DRIFT) && !defined(JUMP)
#ifdef PARA_RES
                #pragma omp parallel for
#endif
                /****** Update residual for all bins (No adaptive method) ******/
                for(j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                        RAND_MESH[j].calculate_second_half(T,DT);             // calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update();
                }
#endif

                sources(RAND_POINTS, DT, N_POINTS);

#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_variables();                   // update the fluid state at vertex
                        RAND_POINTS[i].reset_du();                             // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values();
                        RAND_POINTS[i].con_to_prim();                          // convert these to their corresponding conserved
                }

                if(TBIN_CURRENT == 0){
                        reset_tbins(T, DT, N_TRIANG, N_POINTS, NEXT_DT, RAND_MESH, RAND_POINTS);
                }

// #if defined(FIXED_BOUNDARY) && (defined(NOH) || defined(DF))
//                 for(j=0;j<N_TRIANG;++j){                                         // loop over all triangles in MESH
//                         RAND_MESH[j].check_boundary();                           // calculate flux through TRIANGLE
//                 }
// #endif
                TBIN_CURRENT = (TBIN_CURRENT + 1) % MAX_TBIN;                     // increment time step bin
                T += DT;                                                         // increment time
                l += 1;    // increment step number

        }

        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);

        auto program_end = std::chrono::system_clock::now();
        std::time_t program_end_time = std::chrono::system_clock::to_time_t(program_end);
        std::cout<<"finish computation at "<<std::ctime(&program_end_time);
        std::chrono::duration<double> elapsed_seconds = program_end-program_start;
        std::cout<<"elapsed time: " << elapsed_seconds.count() << "s\n";


        return 0;
}

void set_shock_sensor(int N_POINTS, std::vector<VERTEX> &RAND_POINTS){
    GLOBAL_PRESSURE_MIN = RAND_POINTS[0].get_pressure();
    GLOBAL_PRESSURE_MAX = RAND_POINTS[0].get_pressure();
    double TOTAL_MASS = 0.0;
    for(int i=0;i<N_POINTS;++i){                                       // loop over all vertices
        double PRESSURE =  RAND_POINTS[i].get_pressure();
        if (PRESSURE > GLOBAL_PRESSURE_MAX){GLOBAL_PRESSURE_MAX = PRESSURE;}
        if (PRESSURE < GLOBAL_PRESSURE_MIN){GLOBAL_PRESSURE_MIN = PRESSURE;}
        double POINT_MASS = RAND_POINTS[i].get_dual()*RAND_POINTS[i].get_mass_density();
        double POINT_VELOCITY_SCALE = sqrt(pow(RAND_POINTS[i].get_x_velocity(),2)+pow(RAND_POINTS[i].get_y_velocity(),2));
        GLOBAL_VELOCITY_SCALE += POINT_MASS*POINT_VELOCITY_SCALE;
        TOTAL_MASS += POINT_MASS;
    }
    GLOBAL_VELOCITY_SCALE /= TOTAL_MASS;

    MESH_RESOLUTION = sqrt(SIDE_LENGTH_X*SIDE_LENGTH_Y/N_POINTS*4.0/M_PI);
}
