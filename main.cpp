const int MUL = 1000000;
int convertFloatInt(float _val) { return _val*MUL; }
float convertIntFloat(int _val) { return _val/(float)MUL; }

#include <CGAL/trace.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <vector>
#include <fstream>
// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

using namespace std;

int main()
{
    FT sm_angle = 20.0;
    FT sm_radius = 30;
    FT sm_distance = 0.375;
    PointList points;
    std::ifstream stream("data/kitten.xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(
                              stream,
                              std::back_inserter(points),
                              CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))
    {
      std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
      return EXIT_FAILURE;
    }


//
////    ifstream fin("data/scan6_zoneB.bin", ios::binary);
//    ifstream fin("data/s6s100000.bin", ios::binary);
//    double stationx,stationy,stationz;
//    unsigned int nb_points, donee_supp;
//
//    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
//    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
//    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);
//
//    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);
//    fin.read(reinterpret_cast<char*>(&donee_supp), sizeof donee_supp);
//
//    Chrono chrono = Chrono();
//
//    double x,y,z;
//    for(unsigned int i=0;i</*1000/**/nb_points/**/;i++)
//    {
//
//        fin.read(reinterpret_cast<char*>(&x), sizeof x);
//        fin.read(reinterpret_cast<char*>(&y), sizeof y);
//        fin.read(reinterpret_cast<char*>(&z), sizeof z);
//
//        // TODO à test
////        L.insert(Point(convertFloatInt(x),convertFloatInt(y),convertFloatInt(z)));
////        points.push_back(Point_with_normal(x,y,z,CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())));
//
//        if(i%1000000==0){
//            cout << "complete : " << i << endl;
//        }
//    }

//
//    std::ifstream stream("data/kitten.xyz");
//    if (!stream ||
//    !CGAL::read_xyz_points_and_normals(
//        stream,
//        std::back_inserter(points),
//        CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))
//    {
//        std::cerr << "Error: cannot read file data/kitten.xyz" << std::endl;
//        return EXIT_FAILURE;
//    }

//    Poisson_reconstruction_function function(points.begin(), points.end(),
//    CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()) );
//
//    if ( ! function.compute_implicit_function() )
//        return EXIT_FAILURE;

    return 0;

//    list<Point> L;
//    L.clear();
//    L.push_front(Point(0,0,0));
//    L.push_front(Point(1,0,0));
//    L.push_front(Point(0,1,0));
//    L.push_front(Point(0,0,1));
//    L.push_front(Point(1,1,0));
//    L.push_front(Point(1,0,1));
//    L.push_front(Point(0,1,1));
//    L.push_front(Point(1,1,1));
//
//    Delaunay T(L.begin(), L.end());
//    chrono.printTime();
//    chrono.restart();
//
//    ofstream myfile;
//    myfile.open ("test.txt");
//    myfile << std::fixed;
//    myfile << T;
//    myfile.close();
//    chrono.printTime();
//
//
//    cout << "Triangulation delaunay end" << endl;
//    cout << "Write ply file start" << endl;
//
//    ofstream fout("data/output.ply",  ios::out);
//    fout << "ply" <<endl;
//    fout << "format binary_little_endian 1.0" << endl; // binary
//    fout << "comment author: Etienne Frank University of applied science HES-SO" << endl;
//    fout << "element vertex "<< T.number_of_vertices() << endl;
//    fout << "property float x" << endl;
//    fout << "property float y" << endl;
//    fout << "property float z" << endl;
//    fout << "element face "<< 0/*T.number_of_cells()/**/ << endl;
//    fout << "property list uchar int32 vertex_index" << endl;
//    fout << "end_header" << endl;
//    fout.close();
//    ofstream fout2("data/output.ply",  ios::out| ios::binary| ios::app);
//    ifstream finTemp("test.txt",ios::out);
//
//    string line;
//    getline(finTemp,line);
//    getline(finTemp,line);
//    float xTemp,yTemp,zTemp;
//    for(unsigned int i=0;i<T.number_of_vertices();i++)
//    {
//        finTemp >> xTemp >> yTemp >> zTemp;
//        fout2.write((char *) &xTemp, sizeof(xTemp));
//        fout2.write((char *) &yTemp, sizeof(yTemp));
//        fout2.write((char *) &zTemp, sizeof(zTemp));
//    }
//    getline(finTemp,line);
//    unsigned char listLength = 4;
//    int32_t idx1,idx2,idx3,idx4;
//    for(unsigned int i=0;i<T.number_of_cells();i++)
//    {
//        finTemp >> idx1 >> idx2 >> idx3 >> idx4;
////        fout2.write((char *) &listLength, sizeof(listLength));
////        fout2.write((char *) &idx1, sizeof(idx1));
////        fout2.write((char *) &idx2, sizeof(idx2));
////        fout2.write((char *) &idx3, sizeof(idx3));
////        fout2.write((char *) &idx4, sizeof(idx4));
//        cout << i << " " << idx1 << " " << idx2 << " " << idx3 << " " << idx4 << endl;
//    }
//
//
//    cout << "number of cells" << T.number_of_cells()<<endl;
//    cout << "number of vertices" << T.number_of_vertices()<<endl;
//
//    fout2.close();

//    cout << "PLY file writed" << endl;

    return 0;
}
