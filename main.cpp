const int MUL = 1000000;
int convertFloatInt(float _val)
{
    return _val*MUL;
}
float convertIntFloat(int _val)
{
    return _val/(float)MUL;
}

#include <CGAL/Eigen_solver_traits.h>
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
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <vector>
#include <fstream>
#include <CGAL/pca_estimate_normals.h>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;
typedef std::pair<Point, Vector> PointVectorPair;

using namespace std;

int main()
{
    // Poisson options
    FT sm_angle = 1.0; // Min triangle angle in degrees.
    FT sm_radius = 1; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.1; // Surface Approximation error w.r.t. point set average spacing.
// Reads the point set file in points[].
// Note: read_xyz_points_and_normals() requires an iterator over points
// + property maps to access each point's position and normal.
// The position property map can be omitted here as we use iterators over Point_3 elements.
    PointList points;

    std::list<PointVectorPair> pointsNormals;
//    std::ifstream stream("data/sphere_20k.xyz");
//    if (!stream ||
//            !CGAL::read_xyz_points(stream,
//                                   std::back_inserter(pointsNormals),
//                                   CGAL::First_of_pair_property_map<PointVectorPair>()))
//    {
//        std::cerr << "Error: cannot read file data/sphere_20k.xyz" << std::endl;
//        return EXIT_FAILURE;
//    }
    ifstream fin("data/scan6_zoneB.bin", ios::binary);
//    ifstream fin("data/s6s100000.bin", ios::binary);
    double stationx,stationy,stationz;
    unsigned int nb_points, donee_supp;

    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);

    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);
    fin.read(reinterpret_cast<char*>(&donee_supp), sizeof donee_supp);

//    Chrono chrono = Chrono();

    double x,y,z;
    for(unsigned int i=0; i</*1000/**/nb_points/**/; i++)
    {

        fin.read(reinterpret_cast<char*>(&x), sizeof x);
        fin.read(reinterpret_cast<char*>(&y), sizeof y);
        fin.read(reinterpret_cast<char*>(&z), sizeof z);

        // TODO à test
//        L.insert(Point(convertFloatInt(x),convertFloatInt(y),convertFloatInt(z)));
        PointVectorPair temp;
        temp.first = Point(x,y,z);
        pointsNormals.push_back(temp/*,CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())*/);

        if(i%1000000==0)
        {
            cout << "complete : " << i << endl;
        }
    }


    // Estimates normals direction.
    // Note: pca_estimate_normals() requires an iterator over points
    // as well as property maps to access each point's position and normal.
    const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
    CGAL::pca_estimate_normals(pointsNormals.begin(), pointsNormals.end(),
                               CGAL::First_of_pair_property_map<PointVectorPair>(),
                               CGAL::Second_of_pair_property_map<PointVectorPair>(),
                               nb_neighbors);
    // Orients normals.
    // Note: mst_orient_normals() requires an iterator over points
    // as well as property maps to access each point's position and normal.
    std::list<PointVectorPair>::iterator unoriented_points_begin =
        CGAL::mst_orient_normals(pointsNormals.begin(), pointsNormals.end(),
                                 CGAL::First_of_pair_property_map<PointVectorPair>(),
                                 CGAL::Second_of_pair_property_map<PointVectorPair>(),
                                 nb_neighbors);
    // Optional: delete points with an unoriented normal
    // if you plan to call a reconstruction algorithm that expects oriented normals.
    pointsNormals.erase(unoriented_points_begin, pointsNormals.end());

    // Creates implicit function from the read points using the default solver.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(pointsNormals.begin(), pointsNormals.end(),
                                             CGAL::First_of_pair_property_map<PointVectorPair>(),
                                             CGAL::Second_of_pair_property_map<PointVectorPair>() );

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() )
        return EXIT_FAILURE;

// Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing(pointsNormals.begin(), pointsNormals.end(),
                                                       CGAL::First_of_pair_property_map<PointVectorPair>(),
                                                       6 /* knn = 1 ring */);
// Gets one point inside the implicit surface
// and computes implicit function bounding sphere radius.
    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());
// Defines the implicit surface: requires defining a
// conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);
// Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle, // Min triangle angle (degrees)
            sm_radius*average_spacing, // Max triangle size
            sm_distance*average_spacing); // Approximation error
// Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3, // reconstructed mesh
                            surface, // implicit surface
                            criteria, // meshing criteria
                            CGAL::Manifold_with_boundary_tag()); // require manifold mesh
    if(tr.number_of_vertices() == 0)
        return EXIT_FAILURE;
// saves reconstructed surface mesh
    std::ofstream out("kitten_poisson-20-30-0.375.off");
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
    out << output_mesh;
    return EXIT_SUCCESS;

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
