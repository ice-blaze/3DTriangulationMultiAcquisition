#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include "Chrono.h"

const int MUL = 1000000;
int convertFloatInt(float _val) { return _val*MUL; }
float convertIntFloat(int _val) { return _val/(float)MUL; }

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef Kernel::Point_3 Point;

using namespace std;

int main()
{
// construction from a list of points :
    vector<Point> L;

    ifstream fin("data/scan6_zoneB.bin", ios::binary);
    double stationx,stationy,stationz;
    unsigned int nb_points, donee_supp;

    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);

    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);
    fin.read(reinterpret_cast<char*>(&donee_supp), sizeof donee_supp);

    Chrono chrono = Chrono();

    double x,y,z;
    for(unsigned int i=0;i</*1000/**/nb_points/**/;i++)
    {

        fin.read(reinterpret_cast<char*>(&x), sizeof x);
        fin.read(reinterpret_cast<char*>(&y), sizeof y);
        fin.read(reinterpret_cast<char*>(&z), sizeof z);

        // TODO à test
//        L.insert(Point(convertFloatInt(x),convertFloatInt(y),convertFloatInt(z)));
        L.push_back(Point(convertFloatInt(x),convertFloatInt(y),convertFloatInt(z)));

        if(i%1000000==0){
            cout << "complete : " << i << endl;
        }
    }

//    Delaunay T(L.begin(), L.end());
    chrono.printTime();

//    T.
//    ofstream oFileT("output",ios::out);
//    oFileT << T;
//    Triangulation T1;
//    ifstream iFileT("output",ios::in);
//    return 0;
//}
//    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
//              p != end; ++p)

//    int temp2 = 0;
//    for (Triangulation::Finite_facets_iterator it = T.finite_facets_begin(), end = T.finite_facets_end(); it!= end ; ++it)
//    for (Triangulation::Finite_cells_iterator it = T.finite_cells_begin(), end = T.finite_cells_end(); it!= end ; ++it)
//    {
//        Triangulation::Facet tem = *it;
//        Triangulation::Triangle triangle = T.triangle(tem);
//        it->vertex(0);
//        cout << it->vertex(0)->point() << " " << it->vertex(1)->point() << " " << it->vertex(2)->point() << endl;
//        cout << triangle << endl;
//        return 0;
//        L.
//        cout << it->vertex(0)->info()<<endl;
//        cout << temp2 << endl;
//        temp2++;
//    }
//    cout << T.number_of_edges() << endl;
    cout << "Triangulation delaunay end" << endl;
    cout << "Write ply file start" << endl;

//    ofstream fout("data/output.ply",  ios::out);
//    fout << "ply" <<endl;
//    fout << "format binary_little_endian 1.0" << endl; // binary
//    fout << "comment author: Etienne Frank University of applied science HES-SO" << endl;
//    fout << "element vertex "<< nb_points << endl;
//    fout << "property float x" << endl;
//    fout << "property float y" << endl;
//    fout << "property float z" << endl;
//    fout << "element face "<< 0/*T.number_of_cells()/**/ << endl;
//    fout << "property list uchar int32 vertex_index" << endl;
//    fout << "end_header" << endl;
//    fout.close();
//    ofstream fout2("data/output.ply",  ios::out| ios::binary| ios::app);
//    float temp;
//    for(unsigned int i=0;i<nb_points;i++)
//    {
//        temp = convertIntFloat(L[i].x());
//        fout2.write((char *) &temp, sizeof(temp));
//        temp = convertIntFloat(L[i].y());
//        fout2.write((char *) &temp, sizeof(temp));
//        temp = convertIntFloat(L[i].z());
//        fout2.write((char *) &temp, sizeof(temp));
//    }
//    fout2.close();

    cout << "PLY file writed" << endl;

    return 0;
}
