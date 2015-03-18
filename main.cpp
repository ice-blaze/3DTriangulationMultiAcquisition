#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include "Node.h"


using namespace std;
//const int K = 3;
//typedef struct NOEUD
//{
//    double cle[K];
//    struct NOEUD *g;
//    struct NOEUD *d;
//} NOEUD;
//
//void inserer(NOEUD **arbre, NOEUD *element, int prof)
//{
//    if (*arbre == NULL)
//        (*arbre) = element;
//    else if ((*arbre)->cle[prof%K] < element->cle[prof%K])
//        inserer(&(*arbre)->d, element, prof+1);
//    else
//        inserer(&(*arbre)->g, element, prof+1);
//}
//
//int rapporter(double cle[K], double min[K], double max[K])
//{
//    int i, dedans = 1;
//    for (i = 0; i < K && dedans; i++)
//        dedans = dedans && min[i] <= cle[i] && max[i] >= cle[i];
//    if (dedans) // Rapporter le point
//    {
//        for (i = 0; i < K; i++)
//            printf("%f ", cle[i]);
//        printf("\n");
//    }
//    return dedans;
//}
//void recherche(NOEUD *arbre, double min[K], double max[K],int prof)
//{
//    if (arbre) // arbre non vide
//    {
//        if (arbre->cle[prof%K] >= min[prof%K])
//            recherche(arbre->g, min, max, prof+1);
//        rapporter(arbre->cle, min, max);
//        if (arbre->cle[prof%K] <= max[prof%K])
//            recherche(arbre->d, min, max, prof+1);
//    }
//}

// perf 5.4
int main()
{
    clock_t t1=clock();
    clock_t t2=clock();

    srand (time(NULL));
    Node* racine = new Node();
    racine->key[0] = racine->key[1] = racine->key[2] = 0;

    const long NB_MAX = 60000000;

    Node* noeuds= new Node[60000000];

    for(long i=0; i<NB_MAX; i++)
    {
// NOEUD* noeud = new NOEUD();
// noeud->cle[0] = rand();
// noeud->cle[1] = rand();
// noeud->cle[2] = rand();
// inserer(&racine,noeud,0);
        noeuds[i].key[0] = rand();
        noeuds[i].key[1] = rand();
        noeuds[i].key[2] = rand();
        if(i%10000000==0)
            cout << "complete : " << i << endl;
    }
    t2=clock();
    cout << (float)t2-(float)t1 << " end" << endl;
    return 0;
}

