#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include "Chrono.h"

const int MUL = 1000000;
int convertFloatInt(float _val) { return _val*MUL; }
float convertIntFloat(int _val) { return _val/(float)MUL; }

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
    for(unsigned int i=0;i<1000/*nb_points/**/;i++)
    {

        fin.read(reinterpret_cast<char*>(&x), sizeof x);
        fin.read(reinterpret_cast<char*>(&y), sizeof y);
        fin.read(reinterpret_cast<char*>(&z), sizeof z);

        L.push_front(Point(convertFloatInt(x),convertFloatInt(y),convertFloatInt(z)));

        if(i%1000000==0){
            cout << "complete : " << i << endl;
        }
    }

    Triangulation T(L.begin(), L.end());
    Triangulation::size_type nb_vertices = T.number_of_vertices();

//    T.
//    ofstream oFileT("output",ios::out);
//    oFileT << T;
//    Triangulation T1;
//    ifstream iFileT("output",ios::in);
//    return 0;
//}
//    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
//              p != end; ++p)
    for (Triangulation::All_cells_iterator it = T.all_cells_begin(), end = T.all_cells_end(); it!= end ; ++it)
    {
        for(unsigned int index=0; index < 4; index++) {
//            Triangulation_data_structure::Facet f = std::pair<Cell_handle,int>(it,index);
//            Vertex_triple vt_aux = this->make_vertex_triple(f);
//            Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],vmap[vt_aux.second]);
//            this->make_canonical(vt);
//            inner_map[vt]= f;
        }
    }
    cout << T.number_of_edges() << endl;
    cout << "Triangulation end" << endl;
    chrono.printTime();
    cout << "Write ply file start" << endl;

//    ofstream fout("data/output.ply",  ios::out);
//    fout << "ply" <<endl;
//    fout << "format binary_little_endian 1.0" << endl; // binary
//    fout << "comment author: Etienne Frank University of applied science HES-SO" << endl;
//    fout << "element vertex "<< nb_points << endl;
//    fout << "property float x" << endl;
//    fout << "property float y" << endl;
//    fout << "property float z" << endl;
//    fout << "element face "<< T.number_of_cells() << endl;
//    fout << "property list uchar int32 vertex_index" << endl;
//    fout << "end_header" << endl;
//    fout.close();
//    ofstream fout2("data/output.ply",  ios::out| ios::binary| ios::app);
//    float temp;
//    for(unsigned int i=0;i<nb_points;i++)
//    {
//        temp = convertIntFloat(noeuds[i].cle[0]);
//        fout2.write((char *) &temp, sizeof(temp));
//        temp = convertIntFloat(noeuds[i].cle[1]);
//        fout2.write((char *) &temp, sizeof(temp));
//        temp = convertIntFloat(noeuds[i].cle[2]);
//        fout2.write((char *) &temp, sizeof(temp));
//    }
//    fout2.close();

    cout << "PLY file writed" << endl;

    return 0;
}
