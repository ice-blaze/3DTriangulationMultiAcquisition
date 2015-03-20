#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Triangulation;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Locate_type Locate_type;
typedef Triangulation::Point Point;

using namespace std;

int main()
{
// construction from a list of points :
    list<Point> L;
    L.push_front(Point(0,0,0));
    L.push_front(Point(1,0,0));
    L.push_front(Point(0,1,0));
    Triangulation T(L.begin(), L.end());
    Triangulation::size_type n = T.number_of_vertices();

    ofstream oFileT("output",ios::out);
    oFileT << T;
    Triangulation T1;
    ifstream iFileT("output",ios::in);
    return 0;
}
