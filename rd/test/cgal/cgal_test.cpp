#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       PDT;
typedef CGAL::Iso_rectangle_2<K>                            IR;
typedef PDT::Face_handle                                    Face_handle;
typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Locate_type                                    Locate_type;
typedef PDT::Point                                          Point;
typedef PDT::Iso_rectangle                                  Iso_rectangle;
typedef PDT::Covering_sheets                                Covering_sheets;
int main(){
  float xmax=2.0,ymax=2.0;
  Iso_rectangle domain(0, 0, xmax, ymax); // The cube for the periodic domain

  // construction from a list of points :
  std::list<Point> L;

  int i,j;
  int nx = 100, ny = 10, count = 5000;
  float x,y;

  // for (int i = 0; i < count; ++i){
  //   x = 1.9*(rand() % 10000)/10000.0 + 0.05;
  //   y = 1.9*(rand() % 10000)/10000.0 + 0.05;
  //   L.push_back(Point(x,y));
  // }

  for(i=0;i<nx;++i){
    x = (xmax-0.01*xmax)*float(i)/float(nx);
    for(j=0;j<ny;++j){
      y = (ymax-0.01*ymax)*float(j)/float(ny);
      L.push_back(Point(x,y));
    }
  }

  PDT T(L.begin(), L.end(), domain); // Put the domain with the constructor
  size_t n = T.number_of_vertices();

  T.convert_to_1_sheeted_covering();

  // PDT::Vertex_iterator vit;

  // for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit){
  //   std::cout << n << "\t" << vit->point() << std::endl;
  // }

  // PDT::Face_iterator fit;

  // for (fit = T.faces_begin(); fit != T.faces_end(); ++fit){
  //   std::cout << fit->vertex(0)->point() << std::endl;
  // }

  // std::cout << T << std::endl;

  std::ofstream oFileT("output.tri", std::ios::out);
  // writing file output;
  oFileT << T;

  return 0;
}