#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include "Node.h"


using namespace std;

int main()
{
    ifstream fin("data/scan6_zoneB.bin", ios::binary);
    double stationx,stationy,stationz;
    unsigned int nb_points, donee_supp;

    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);

    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);
    fin.read(reinterpret_cast<char*>(&donee_supp), sizeof donee_supp);
    clock_t t1=clock();
    clock_t t2=clock();

    srand (time(NULL));
    Node* racine = new Node();
    racine->key[0] = racine->key[1] = racine->key[2] = 0;

//    const long NB_MAX = 60000000;

    Node* noeuds= new Node[nb_points];

    for(long i=0; i<nb_points; i++)
    {
        fin.read(reinterpret_cast<char*>(&noeuds[i].key[0]), sizeof noeuds[i].key[0]);
        fin.read(reinterpret_cast<char*>(&noeuds[i].key[1]), sizeof noeuds[i].key[1]);
        fin.read(reinterpret_cast<char*>(&noeuds[i].key[2]), sizeof noeuds[i].key[2]);

        Node::add(&racine,&noeuds[i]);

        if(i%10000000==0)
            cout << "complete : " << i << endl;
    }
    t2=clock();
    cout << (float)t2-(float)t1 << " end" << endl;
    return 0;
}

