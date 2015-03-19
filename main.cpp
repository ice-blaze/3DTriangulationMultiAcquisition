#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cstdlib>

using namespace std;

const int K = 3;

typedef struct NOEUD
{
    double cle[K];
    struct NOEUD *g;
    struct NOEUD *d;
} NOEUD;

void inserer(NOEUD **arbre, NOEUD *element, int prof)
{
    if (*arbre == NULL)
        (*arbre) = element;
    else if ((*arbre)->cle[prof%K] < element->cle[prof%K])
        inserer(&(*arbre)->d, element, prof+1);
    else
        inserer(&(*arbre)->g, element, prof+1);
}

int rapporter(double cle[K], double min[K], double max[K])
{
    int i, dedans = 1;
    for (i = 0; i < K && dedans; i++)
        dedans = dedans && min[i] <= cle[i] && max[i] >= cle[i];
    if (dedans) // Rapporter le point
    {
        for (i = 0; i < K; i++)
            printf("%f ", cle[i]);
        printf("\n");
    }
    return dedans;
}

void recherche(NOEUD *arbre, double min[K], double max[K],int prof)
{
    if (arbre) // arbre non vide
    {
        if (arbre->cle[prof%K] >= min[prof%K])
            recherche(arbre->g, min, max, prof+1);
        rapporter(arbre->cle, min, max);
        if (arbre->cle[prof%K] <= max[prof%K])
            recherche(arbre->d, min, max, prof+1);
    }
}

int main()
{
    srand (time(NULL));

    NOEUD* racine = new NOEUD();
    racine->cle[0] = racine->cle[1] = racine->cle[2] = 0;

    std::ifstream fin("data/scan6_zoneB.bin", ios::binary);
    double stationx,stationy,stationz;
    unsigned int nb_points, donee_supp;

    fin.read(reinterpret_cast<char*>(&stationx), sizeof stationx);
    fin.read(reinterpret_cast<char*>(&stationy), sizeof stationy);
    fin.read(reinterpret_cast<char*>(&stationz), sizeof stationz);

    fin.read(reinterpret_cast<char*>(&nb_points), sizeof nb_points);
    fin.read(reinterpret_cast<char*>(&donee_supp), sizeof donee_supp);

    clock_t t1=clock();
    NOEUD* noeuds= new NOEUD[nb_points];

    int x,y,z;
    for(unsigned int i=0;i<nb_points;i++)
    {
        fin.read(reinterpret_cast<char*>(&noeuds[i].cle[0]), sizeof noeuds[i].cle[0]);
        fin.read(reinterpret_cast<char*>(&noeuds[i].cle[1]), sizeof noeuds[i].cle[1]);
        fin.read(reinterpret_cast<char*>(&noeuds[i].cle[2]), sizeof noeuds[i].cle[2]);

        inserer(&racine,&noeuds[i],0);

        if(i%10000000==0)
            cout << "complete : " << i << endl;
    }
    clock_t t2=clock();
    cout << (float)t2-(float)t1 << " end" << endl;
    return 0;
}
