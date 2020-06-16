#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>       Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>             P3DT3;
typedef P3DT3::Point             Point;
typedef P3DT3::Iso_cuboid        Iso_cuboid;
typedef P3DT3::Vertex_handle     Vertex_handle;
typedef P3DT3::Cell_handle       Cell_handle;
typedef P3DT3::Locate_type       Locate_type;


// #define RANDOMIC
#define UNIFORMIC
// #define UNIFORMOFFSETIC
// #define PERTUNIFORMOFFSETIC

int main(){
        float xmax=1.0,ymax=1.0,zmax=1.0;
        Iso_cuboid domain(0, 0, 0, xmax, ymax, zmax); // The cube for the periodic domain

        // construction from a list of points :
        std::list<Point> L;

        int i,j,k;
        int nx=4, ny=4, nz=4, count=nx*ny*nz;
        float x,y,z,xmove,ymove,zmove;

#ifdef RANDOMIC
        for (int i = 0; i < count; ++i){
                x = xmax*(rand() % 10000)/10000.0;
                y = ymax*(rand() % 10000)/10000.0;
                z = zmax*(rand() % 10000)/10000.0;
                // std::cout << x << "\t" << y << "\t" << z << std::endl;
                L.push_back(Point(x,y,z));
        }
#endif

#ifdef UNIFORMIC
        for(i=0; i < (nx); ++i){
                x = (xmax) * float(i) / float(nx);
                for(j=0; j < (ny); ++j){
                        y = (ymax) * float(j) / float(ny);
                        for(k=0; k < (nz); ++k){
                                z = (zmax) * float(k) / float(nz);
                                L.push_back(Point(x,y,z));
                        }
                }
        }
#endif

        P3DT3 T(L.begin(), L.end(), domain); // Put the domain with the constructor

        T.convert_to_1_sheeted_covering();

        std::ofstream oFileT("output.txt", std::ios::out);
        // writing file output;
        oFileT << T;

        return 0;
}
