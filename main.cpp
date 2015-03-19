#include <iostream>
#include <fstream>

using namespace std;



int main()
{
    std::ifstream fin("data/scan6_zoneB.bin", ios::binary);
    double stationx,stationy,stationz;
    unsigned int nb_points, donee_supp;

    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);

    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);

//    fin >> stationx;
//    fin >> stationy;
//    fin >> stationz;

//    fin >> nb_points;
//    fin >> donee_supp;
//    std::cout << d << std::endl;
    cout << "station : " << stationx << " " << stationy << " " << stationz << endl;
    cout << "nb points " << nb_points << endl;
    return 0;
}
