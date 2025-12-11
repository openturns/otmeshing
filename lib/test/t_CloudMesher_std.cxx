
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

int main(int argc, char* argv[])
{
  double pointsIn[][4] = {
    { 42.89, 0, 60.55, 30.72 },
    { 45.65, 50.83, 50.37, 16.13 },
    { 79.06, 57.84, 61.59, 2.52 },
    { 44.47, 39.46, 39.53, 28.72},
    { 0, 100, 0, 0 },
    { 66.95, 100, 33.6, 0 },
    { 42.89, 0, 0, 30.72 },
    { 100, 100, 100, 100 }
  };
 
  // typedef CGAL::Delaunay_triangulation<CGAL::Epick_d< CGAL::Dimension_tag<7> > >      T;
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d< CGAL::Dynamic_dimension_tag > >      T;
  
  const int n = 3 + argc;
  std::cout << "n=" << n << std::endl;
  T dt(n);

  std::vector<T::Point> points;
  points.reserve(8);
  for (int j = 0; j < 8; ++j) {
    T::Point p(&pointsIn[j][0], &pointsIn[j][4]);
    points.push_back(p);
  }
 
  // T::Vertex_handle hint;
  // int i = 0;
  for (std::vector<T::Point>::iterator it = points.begin(); it != points.end(); ++it)
    dt.insert(*it);
    // }
    // printf("Processing: %d/%d\n", ++i, (int)points.size());
  // }

  printf("number_of_vertices=%ld\n", dt.number_of_vertices());

  // UnsignedInte
  for (T::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); ++ vi)
    {
      std::cout << "eee=" << vi->point() << std::endl;
      // vertices.add(Point(vi->point().cartesian_begin(), vi->point().cartesian_end()));
      // vertexToIndexMap[vi->point()] = iV;
      // ++ iV;
    }
  
  return 0;
}

