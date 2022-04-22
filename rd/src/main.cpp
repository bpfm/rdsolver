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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       PDT;
typedef CGAL::Iso_rectangle_2<K>                            IR;
//typedef PDT::Face_handle                                    Face_handle;
typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Locate_type                                    Locate_type;
typedef PDT::Point                                          Point;
typedef PDT::Iso_rectangle                                  Iso_rectangle;
typedef PDT::Covering_sheets                                Covering_sheets;

typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT>    Vbb;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, GT, Vbb>  Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2<GT>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                    Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds>              Delaunay;
typedef Delaunay::Point                                                 Point;
typedef Delaunay::Face_handle                                       Face_handle;



#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#endif
#ifdef SEDOV
    int POINT_CHECK = 0;
#endif



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
        int SNAP_MESH_ID = 0;
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

        //create output directory
        fs::path sourceFile = "src/constants.h";
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



#ifdef GENERATE_MESH   //(not used)
//    float xmax=1.0,ymax=1.0;
//    Iso_rectangle domain(0, 0, xmax, ymax); // The cube for the periodic domain
//
//    // construction from a list of points :
//    std::list<std::pair<Point,unsigned>> PointList;
//    //std::list<Point> PointList;
//
//    int nx=40, ny=40, count=nx*ny;
//    float x,y,xmove,ymove;
//
//    #ifdef UNIFORMIC
//        unsigned int point_count = 0;
//            for(i=0; i < (nx); ++i){
//                    x = (xmax) * float(i) / float(nx);
//                    for(j=0; j < (ny); ++j){
//                            y = (ymax) * float(j) / float(ny);
//                            //L.push_back(Point(x,y));
//                            PointList.push_back(std::make_pair( Point(x,y), point_count));
//                            point_count += 1;
//                    }
//            }
//            cout<<"point_count: "<<point_count<<endl;
//    #endif
//
//    //PDT T; // Put the domain with the constructor
//    Delaunay Delaunay_Triangulation;
//    Delaunay_Triangulation.insert(PointList.begin(), PointList.end());
//    size_t n_vertex = Delaunay_Triangulation.number_of_vertices();
//    size_t n_triangle = Delaunay_Triangulation.number_of_faces();
//    Delaunay_Triangulation.convert_to_1_sheeted_covering();
//
//    std::ofstream oFileT("./test.txt", std::ios::out);
//    //writing file output;
//    oFileT << Delaunay_Triangulation;

#endif

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
                        printf("%d\t%f\t%f\t%f\n", POINT_CHECK, PRESSURE_AIM, RAND_POINTS[i].get_pressure(), ETOT);
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
        NEXT_DT = 0.0;                                                            // set first timestep to zero

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
                        write_active(RAND_MESH, N_TRIANG, SNAP_ID, TBIN_CURRENT);
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        SNAP_ID ++;

                #ifdef MOVING_MESH
                        write_mesh(RAND_POINTS,RAND_MESH,T,N_POINTS,N_TRIANG,SNAP_MESH_ID);
                        SNAP_MESH_ID ++;
                #endif
                }



        /****** 1st order update ***************************************************************************************************/

                //dual_halfdt is used to update resiual; while dual is used for CFL;
                //before updating residual, we reset dual_halfdt and then calculate the new one.
#ifdef PARA_UP
#pragma omp parallel for
#endif
                for(i=0;i<N_POINTS;++i){
                        RAND_POINTS[i].set_dual_halfdt(0.0);
                }


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
                /****** setup normal vectors for n+1/2 timestep, and calculate n+1/2 dual areas for vertices ******/
                /****** warning: Before calculating residuals, we must set up normals for all triangles. ******/
                for(j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        RAND_MESH[j].setup_normals(DT);
                }


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

                //moving mesh: update positions of vertices
#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for (i=0;i<N_POINTS;++i){
                        double NEW_X = RAND_POINTS[i].get_x() + RAND_POINTS[i].get_sigma_x()*DT;
                        NEW_X = std::fmod(NEW_X,SIDE_LENGTH_X);
                        RAND_POINTS[i].set_x(NEW_X);
                        double NEW_Y = RAND_POINTS[i].get_y() + RAND_POINTS[i].get_sigma_y()*DT;
                        NEW_Y = std::fmod(NEW_Y,SIDE_LENGTH_Y);
                        RAND_POINTS[i].set_y(NEW_Y);


                        RAND_POINTS[i].set_dual(0.0);
                        //setup new mesh velocity
                        RAND_POINTS[i].set_sigma_x(1.5);
                        RAND_POINTS[i].set_sigma_y(0.0);
                }
#ifdef PARA_UP
                #pragma omp parallel for
#endif
                for (j=0;j<N_TRIANG;++j){
                        RAND_MESH[j].setup_positions_area();
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
