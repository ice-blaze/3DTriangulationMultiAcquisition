/***********************************************************************
But: Calcul d'une triangulation de Delaunay sur des mesures effectuées
par une station laser dont on connaît l'emplacement, génération d'un
fichier au format ply

Fichier: triangule_scan_ply.c
Auteur & copyright: É. Taillard, 17.4.2015
Compilation : gcc -Wall -O3 triangule_scan_ply.c -lm
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <inttypes.h> // coordonnées des points entières de type int32_t
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <chrono>
#include "sys/types.h"
#include "sys/sysinfo.h"

// Pour passer de coordonnées de type double en int32_t
#define PRECISION 10000000

#include <bitset>
#include <map>
#include <algorithm>
#include <iostream>
using namespace std;

const unsigned char nb_cote_triangle = 3; //Pour l'écriture de 3 en binaire

// Les points de mesure à prendre en considération seront à l'intérieur d'un
// parallélipipède compris entre [(xmin, ymin, zmin), (xmax, ymax, zmax)]
//double xmin = -1.0, xmax = 3.0, ymin = -2.0, ymax = 3.0, zmin = -1.0, zmax = 3.0;
//double xmin = -10.0, xmax = 30.0,
//       ymin = -20.0, ymax = 30.0,
//       zmin = -10.0, zmax = 30.0;

const double MAX = 10000.0;
double xmin = -MAX, ymin = -MAX, zmin = -MAX,
       xmax = MAX,  ymax = MAX,  zmax = MAX;

bool is_equal(double d1, double d2){
  int i1 = d1*PRECISION;
  int i2 = d2*PRECISION;
  if(i1==i2) return true;
  return false;
}

struct POINT {double x; double y; double z;unsigned int stations;
  POINT() {stations=0;}
  POINT(double x, double y, double z) : x(x), y(y), z(z) {stations=0;}
  bool operator<(const POINT &o) const {
    if (!is_equal(x,o.x)) { return x < o.x; }
    if (!is_equal(y,o.y)) { return y < o.y; }
    return z < o.z;
  }
  bool operator==(const POINT &o) const {
    if (is_equal(x,o.x) && is_equal(y,o.y) && is_equal(z,o.z))
      return true;
    return false;
  }
  bool operator!=(const POINT &o) const {
    if (is_equal(x,o.x) && is_equal(y,o.y) && is_equal(z,o.z))
      return false;
    return true;
  }
  friend POINT operator+(POINT lhs, const POINT &o) {
    lhs.x += o.x; lhs.y += o.y; lhs.z += o.z;
    return lhs;
  }
  friend POINT operator/(POINT lhs, const double scalaire) {
    lhs.x /= scalaire; lhs.y /= scalaire; lhs.z /= scalaire;
    return lhs;
  }
};
/* DEBUG */
void print(POINT p) {cout << "point : " << p.x << " " << p.y << " " << p.z << endl;}
typedef struct {unsigned char x; unsigned char y; unsigned char z;} COULEUR;

// puissance de 2
int pow2( unsigned int pwr ) {
   return 1 << pwr;
}

// est une puissance de 2 ?
bool isPow2 (unsigned int x){
  return ((x != 0) && ((x & (~x + 1)) == x));
}

const int FALSE = 0;
const int TRUE = 1;

// L'arborescence a environ 9 triangles par point; 10 pour être confortable
const int nb_triangles_par_point =10;
int nb_max_triangles;

double seconds(clock_t cpu)
{ return (double)(clock() - cpu) / CLOCKS_PER_SEC;}

typedef struct{
  int p1, p2, p3,           // Points du triangle
      t1, t2, t3,           // Triangles adjacents ou triangles fils
      final;}  t_triangle;  // Si le triangle n'est pas divisé


/*************** L'Ecuyer random number generator ***************/
double rando(void)
 {
  static int x10 = 12345, x11 = 67890, x12 = 13579, /* initial value*/
             x20 = 24680, x21 = 98765, x22 = 43210; /* of seeds*/
  const int m = 2147483647; const int m2 = 2145483479;
  const int a12= 63308; const int q12=33921; const int r12=12979;
  const int a13=-183326; const int q13=11714; const int r13=2883;
  const int a21= 86098; const int q21=24919; const int r21= 7417;
  const int a23=-539608; const int q23= 3976; const int r23=2071;
  const double invm = 4.656612873077393e-10;
  int h, p12, p13, p21, p23;
  h = x10/q13; p13 = -a13*(x10-h*q13)-h*r13;
  h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
  if (p13 < 0) p13 = p13 + m; if (p12 < 0) p12 = p12 + m;
  x10 = x11; x11 = x12; x12 = p12-p13; if (x12 < 0) x12 = x12 + m;
  h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
  h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
  if (p23 < 0) p23 = p23 + m2; if (p21 < 0) p21 = p21 + m2;
  x20 = x21; x21 = x22; x22 = p21-p23; if(x22 < 0) x22 = x22 + m2;
  if (x12 < x22) h = x12 - x22 + m; else h = x12 - x22;
  if (h == 0) return(1.0); else return(h*invm);
 }

/*********** return an integer between low and high *************/
int unif(int low, int high)
 {return low + (int)((double)(high - low + 1) * rando()) ;}

void transpose(int *a, int *b) {int temp = *a; *a = *b; *b = temp;}

/***************************** Delaunay **************************/

//Indique si le point (px, py) est à gauche du segment p1->p2
int a_gauche(int32_t px, int32_t py, int p1, int p2, int32_t* x, int32_t* y)
{ int32_t wx, wy, vx, vy;
  int64_t vxwy, vywx, vxwx, vywy;
  wx = x[p2] - x[p1];
  wy = y[p2] - y[p1];
  vx = px - x[p1];
  vy = py - y[p1];
  vxwy = (int64_t)vx * (int64_t)wy;
  vywx = (int64_t)vy * (int64_t)wx;
  vxwx = (int64_t)vx * (int64_t)wx;
  vywy = (int64_t)vy * (int64_t)wy;
  if (vxwy > vywx)
    return -1;  // à droite
  else if ((vxwy == vywx) && (vxwx >= 0) && (vywy >= 0) &&
           (abs (vx) <= abs(wx)) && (abs (vy) <= abs(wy)) )
    return 0;  // sur le segment
  else
    return 1; // à gauche
}

//Indique si le point (px, py) est dans un triangle donné
int dans_triangle(int32_t px, int32_t py, t_triangle triangle,
                  int32_t* x, int32_t* y)
{
  int est_a_gauche;

  est_a_gauche = a_gauche(px, py, triangle.p1, triangle.p2, x, y);
  if (est_a_gauche < 0) return FALSE;
  else if (est_a_gauche == 0) return triangle.t1;

  est_a_gauche = a_gauche(px, py, triangle.p2, triangle.p3, x, y);
  if (est_a_gauche < 0) return FALSE;
  else if (est_a_gauche == 0) return triangle.t2;

  est_a_gauche = a_gauche(px, py, triangle.p3, triangle.p1, x, y);
  if (est_a_gauche < 0) return FALSE;
  else if (est_a_gauche == 0) return triangle.t3;

  return TRUE;
}

//Identifie dans quel triangle le point (px, py) se trouve
//Si le point est sur une arête, il est à la fois dans dedans1 et dedans2
void identifie_triangle(int32_t px, int32_t py, t_triangle * arbre,
                        int32_t* x, int32_t* y, int * dedans1, int * dedans2)
{ int triangle = 0, tri2;
  if (!dans_triangle(px, py, arbre[triangle], x, y) )
       printf("Point pas dans le triangle d'origine, %d %d\n", px, py);
  *dedans1 = 0; *dedans2 = -1;
  while (!arbre[triangle].final)
  { if ( (tri2 = dans_triangle(px, py, arbre[arbre[triangle].t1], x, y)) )
    { //printf("dans premier fils\n");
      if (arbre[arbre[triangle].t1].final)
      { *dedans1 = arbre[triangle].t1;
        if (tri2 > 1) *dedans2 = tri2;
      }
      triangle = arbre[triangle].t1;
    }

    else if ( (tri2 = dans_triangle(px, py, arbre[arbre[triangle].t2], x, y)) )
    { //printf("dans second fils\n");
      if (arbre[arbre[triangle].t2].final)
      { *dedans1 = arbre[triangle].t2;
        if (tri2 > 1) *dedans2 = tri2;
      }
      triangle = arbre[triangle].t2;
    }

    else if (arbre[triangle].t3 > 0 &&
             (tri2 = dans_triangle(px, py, arbre[arbre[triangle].t3], x, y)) )
    { //printf("dans dernier fils\n");
      if (arbre[arbre[triangle].t3].final)
      { *dedans1 = arbre[triangle].t3;
        if (tri2 > 1) *dedans2 = tri2;
      }
      triangle = arbre[triangle].t3;
    }

    else
    {
      printf("(%d, %d) Dans le triangle %d (%d) mais pas dans ses fils %d %d %d\n",
               px, py, triangle, dans_triangle(px, py, arbre[triangle], x, y),
               arbre[triangle].t1, arbre[triangle].t2, arbre[triangle].t3);
      return;
    }
  }
  if (*dedans2 > 0)
  {
    triangle = *dedans2;
    while (!arbre[triangle].final)
    {if (dans_triangle(px, py, arbre[arbre[triangle].t1], x, y) )
       triangle = arbre[triangle].t1;
     else if (dans_triangle(px, py, arbre[arbre[triangle].t2], x, y) )
       triangle = arbre[triangle].t2;
     else if (dans_triangle(px, py, arbre[arbre[triangle].t3], x, y) )
       triangle = arbre[triangle].t3;
     else
       printf("dans %d mais pas dans ses fils\n", triangle);
    }
    *dedans2 = triangle;
  }
}

double sqr (double x) {return x * x;}

/* vérifie si le point (x4, y4) est en dehors du cercle circonscrit par le
   triangle 1-2-3; on fait l'hypothèse que les imprécisions des doubles
   sont acceptables lorqu'on considère une arête comme légale */
int legal(int32_t x1, int32_t y1, int32_t x2, int32_t y2,
          int32_t x3, int32_t y3, int32_t x4, int32_t y4)
{ double x5, y5, r; // centre du cercle, rayon au carré
  double ax, ay, bx, by, cx, cy, dx, dy; // perpendiculaires aux segments 1-2 et 1-3
  double denominateur;  // dénominateur pour le calcul de l'intersection des segments

  ax = (x1 + x2)/2.0; ay = (y1 + y2)/2.0;
  bx = ax -  (y2 - y1); by = ay +  (x2 - x1);
  cx = (x1 + x3)/2.0; cy = (y1 + y3)/2.0;
  dx = cx - (y3 - y1); dy = cy +  (x3 - x1);
  denominateur =  (ax - bx) * (cy - dy) - (ay - by)*(cx - dx);
  x5 = (ax * by - ay * bx) * (cx - dx) - (ax - bx) * (cx * dy - cy * dx);
  x5 /= denominateur;
  y5 = (ax * by - ay * bx) * (cy - dy) - (ay - by) * (cx * dy - cy * dx);
  y5 /= denominateur;
  r = sqr(x1 - x5) + sqr(y1 - y5);

  if ( r > sqr(x4 - x5) + sqr(y4 - y5) )
  {
    return FALSE;
  }
  else
    return TRUE;
}

//Découpe le triangle t en trois triangles t1, t2 et t3
void decoupe_triangle(int t, int t1, int t2, int t3, t_triangle* arbre)
{
  arbre[t].t1 = t1; arbre[t].t2 = t2; arbre[t].t3 = t3; arbre[t].final = FALSE;
}

//Modifie la numérotation du triangle adjacent après une découpe
//Le triangle qui était adjacent à "avant" devient adjacent à "apres"
void modifie_triangle_adjacent(int adjacent, int avant, int apres,
                               t_triangle* arbre)
{  if (adjacent > 0)
  { if      (arbre[adjacent].t1 == avant) arbre[adjacent].t1 = apres;
    else if (arbre[adjacent].t2 == avant) arbre[adjacent].t2 = apres;
    else if (arbre[adjacent].t3 == avant) arbre[adjacent].t3 = apres;
    else {printf("%d pas dans les triangles adjacents de %d\n", avant, adjacent);
    return;}
  }
}


//Légalise l'arête commune aux triangles t1 et t2
void legaliser(int t1, int t2, int32_t* x, int32_t* y,
               int* nb_triangles, t_triangle* arbre)
{
 int pt1, pt2, pt3, pt4, t3, t4, t5, t6, t7, t8;
 if (legal(x[arbre[t2].p1], y[arbre[t2].p1],
           x[arbre[t2].p2], y[arbre[t2].p2],
           x[arbre[t2].p3], y[arbre[t2].p3],
           x[arbre[t1].p1], y[arbre[t1].p1]) )
   return; // arête légale, on ne fait rien
 else
 {if (arbre[t1].t2 != t2)
   {printf("t2 pas egal au second triangle adjacent de t1\n");
    return;
   }
  pt1 = arbre[t1].p1;
  pt2 = arbre[t1].p2;
  pt3 = arbre[t1].p3;
  t3 = arbre[t1].t1;
  t4 = arbre[t1].t3;
  if (arbre[t2].t1 == t1)
  { pt4 = arbre[t2].p3;
    t5 = arbre[t2].t2; t6 = arbre[t2].t3;
  }
  else if (arbre[t2].t2 == t1)
  { pt4 = arbre[t2].p1;
    t5 = arbre[t2].t3; t6 = arbre[t2].t1;
  }
  else if (arbre[t2].t3 == t1)
  { pt4 = arbre[t2].p2;
    t5 = arbre[t2].t1; t6 = arbre[t2].t2;
  }
  else
  { printf("t1 pas dans triangles adjacents de t2\n");
    return;
  }
  t7 = *nb_triangles;
  t8 = *nb_triangles + 1;

  *nb_triangles = *nb_triangles + 2;

  decoupe_triangle(t1, t7, t8, -1, arbre);
  arbre[t7] = (t_triangle){pt1, pt2, pt4, t3, t5, t8, TRUE};

  decoupe_triangle(t2, t7, t8, -1, arbre);
  arbre[t8] = (t_triangle){pt1, pt4, pt3, t7, t6, t4, TRUE};

  modifie_triangle_adjacent(t3, t1, t7, arbre);
  modifie_triangle_adjacent(t4, t1, t8, arbre);
  modifie_triangle_adjacent(t5, t2, t7, arbre);
  modifie_triangle_adjacent(t6, t2, t8, arbre);

  if (t5 > 0) legaliser(t7, t5, x, y, nb_triangles, arbre);
  if (t6 > 0) legaliser(t8, t6, x, y, nb_triangles, arbre);
 }
}

//Insère le point d'indice p dans la triangulation
int insere_point(int p, int32_t*x, int32_t* y, int* nb_triangles,
                 t_triangle* arbre)
{int dedans1, dedans2;
 int pt1, pt2, pt3, pt4,
     t1, t2, t3, t4, t5, t6, t7, t8;
 identifie_triangle(x[p], y[p], arbre, x, y, & dedans1,  & dedans2);

 if (dedans1 >= 0 && dedans2 >= 0) // point superposé ou sur un arête
 { if ( (x[arbre[dedans1].p1] == x[p] && y[arbre[dedans1].p1] == y[p]) ||
        (x[arbre[dedans1].p2] == x[p] && y[arbre[dedans1].p2] == y[p]) ||
        (x[arbre[dedans1].p3] == x[p] && y[arbre[dedans1].p3] == y[p]) )
    // Le point est  ignoré car superposé à un point déjà existant
    return 0;
    else // point sur une arête mais pas superposé
    {
       if (arbre[dedans1].t1 == dedans2)
       {
         pt1 = arbre[dedans1].p1;
         pt2 = arbre[dedans1].p2;
         pt3 = arbre[dedans1].p3;
         t2 = arbre[dedans1].t2; t3 = arbre[dedans1].t3;
       }
       else if (arbre[dedans1].t2 == dedans2)
       {
         pt1 = arbre[dedans1].p2;
         pt2 = arbre[dedans1].p3;
         pt3 = arbre[dedans1].p1;
         t2 = arbre[dedans1].t3; t3 = arbre[dedans1].t1;
       }
       else if (arbre[dedans1].t3 == dedans2)
       {
         pt1 = arbre[dedans1].p3;
         pt2 = arbre[dedans1].p1;
         pt3 = arbre[dedans1].p2;
         t2 = arbre[dedans1].t1; t3 = arbre[dedans1].t2;
       }
       else
       { printf("dedans2 pas dans les triangle adjacents de dedans1\n");
         return 0;
       }

       if (arbre[dedans2].t1 == dedans1)
       {
         pt4 = arbre[dedans2].p3;
         t1 = arbre[dedans2].t2; t4 = arbre[dedans2].t3;
       }
       else if (arbre[dedans2].t2 == dedans1)
       {
         pt4 = arbre[dedans2].p1;
         t1 = arbre[dedans2].t3; t4 = arbre[dedans2].t1;
       }
       else if (arbre[dedans2].t3 == dedans1)
       {
         pt4 = arbre[dedans2].p2;
         t1 = arbre[dedans2].t1; t4 = arbre[dedans2].t2;
       }
       else
       { printf("dedans1 pas dans les triangle adjacents de dedans2\n");
         return 0;
       }

       t5 = *nb_triangles;
       t6 = *nb_triangles + 1;
       t7 = *nb_triangles + 2;
       t8 = *nb_triangles + 3;

       (*nb_triangles) += 4;

       decoupe_triangle(dedans1, t5, t6, -1, arbre);
       arbre[t5] = (t_triangle){p, pt2, pt3, t8, t2, t6, TRUE};
       modifie_triangle_adjacent(t2, dedans1, t5, arbre);
       arbre[t6] = (t_triangle){p, pt3, pt1, t5, t3, t7, TRUE};
       modifie_triangle_adjacent(t3, dedans1, t6, arbre);

       decoupe_triangle(dedans2, t7, t8, -1, arbre);
       arbre[t7] = (t_triangle){p, pt1, pt4, t6, t1, t8, TRUE};
       modifie_triangle_adjacent(t1, dedans2, t7, arbre);
       arbre[t8] = (t_triangle){p, pt4, pt2, t7, t4, t5, TRUE};
       modifie_triangle_adjacent(t4, dedans2, t8, arbre);

       if (t2 > 0) legaliser(t5, t2, x, y, nb_triangles, arbre);
       if (t3 > 0) legaliser(t6, t3, x, y, nb_triangles, arbre);
       if (t1 > 0) legaliser(t7, t1, x, y, nb_triangles, arbre);
       if (t4 > 0) legaliser(t8, t4, x, y, nb_triangles, arbre);

    } // point sur une arête mais pas superposé
  } //if (dedans1 >= 0 && dedans2 >= 0)

  else if (dedans1 >= 0) // point à l'intérieur du triangle dedans1
  {
     pt1 = arbre[dedans1].p1;
     pt2 = arbre[dedans1].p2;
     pt3 = arbre[dedans1].p3;
     t1 = arbre[dedans1].t1;
     t2 = arbre[dedans1].t2;
     t3 = arbre[dedans1].t3;
     t4 = *nb_triangles;
     t5 = *nb_triangles + 1;
     t6 = *nb_triangles + 2;

     (*nb_triangles) += 3;

     decoupe_triangle(dedans1, t4, t5, t6, arbre);
     arbre[t4] = (t_triangle){p, pt1, pt2, t6, t1, t5, TRUE};
     modifie_triangle_adjacent(t1, dedans1, t4, arbre);
     arbre[t5] = (t_triangle){p, pt2, pt3, t4, t2, t6, TRUE};
     modifie_triangle_adjacent(t2, dedans1, t5, arbre);
     arbre[t6] = (t_triangle){p, pt3, pt1, t5, t3, t4, TRUE};
     modifie_triangle_adjacent(t3, dedans1, t6, arbre);

     if (t1 > 0) legaliser(t4, t1, x, y, nb_triangles, arbre);
     if (t2 > 0) legaliser(t5, t2, x, y, nb_triangles, arbre);
     if (t3 > 0) legaliser(t6, t3, x, y, nb_triangles, arbre);

  }
  else
  {printf("Triangle non identifie\n");
   return 0;
  }
  return 1;
}

//Construit une triangulation de Delaunay sur un semi de n points
//indicés de 0 à n-1. Les points n, n+1 et n+2 doivent englober la scène.
//Retourne le nombre de triangles totaux de la triangulation, y compris
//les triangles intermédiaires de la structure de données permettant
//d'identifier rapidement le triangle de la triangulation dans lequel un
//point se trouve
int construit_delaunay(int n, int32_t* x, int32_t* y, t_triangle* arbre)
{
  int* permutation;
  int i, nb_triangles;
  permutation = (int*) malloc( n*sizeof(int));
  for (i = 0; i < n; i++) permutation[i] = i;
  for (i = 0; i < n-1; i++)
    transpose(&permutation[i], &permutation[unif(i, n-1)]);

  arbre[0] =(t_triangle){n, n+1, n+2, -1, -1, -1, TRUE};
  nb_triangles = 1;

  int sum =0;
  for (i = 0; i < n; i++)
  {//TODO insere point bug avec la table
     if (!insere_point(permutation[i], x, y, &nb_triangles, arbre))
//       printf("%de point pas insere; %d;\n", i, nb_triangles);
        ++sum;
     if (nb_triangles >= nb_triangles_par_point*n)
     {
       printf("Arbre mal dimensionne; augmenter nb_triangles_par_point\n");
       return 0;
     }
     if (i%(n/100) == 0)
     { printf("\b\b\b%2d%%", 100*i/n); fflush(stdout);}
  }

  printf(" nb de points ignore´ %d\n",sum);

  free(permutation);
  return nb_triangles;
}

/***************************** fin delaunay ************************/

//Calcule l'azimuth et l'élévation d'un point p vu de "origine"
void calcule_phi_theta(POINT origine, POINT p, int32_t *phi, int32_t * theta)
{
  double proj_hor = sqrt( (p.x - origine.x)*(p.x - origine.x) +
                          (p.y - origine.y)*(p.y - origine.y) );
  *phi = PRECISION * atan2(p.y - origine.y, p.x - origine.x) + 0.5;
  *theta = PRECISION * atan2(p.z - origine.z, proj_hor) + 0.5;
}

// recoit la position de la station du fichier
POINT extrait_station(const char chemin_ply[]){

  FILE* fichier = fopen(chemin_ply, "rb");
  POINT origine;
  fread(&origine, sizeof origine, 1, fichier);
  fclose(fichier);
  return origine;
}

/* le point p est du même côté du triangle p1 p2 p3 que le point q */
int visible(POINT p, POINT q, POINT p1, POINT p2, POINT p3)
{ POINT v12 = {p2.x-p1.x, p2.y-p1.y, p2.z-p1.z},
        v13 = {p3.x-p1.x, p3.y-p1.y, p3.z-p1.z},
        v1p = {p.x-p1.x, p.y-p1.y, p.z-p1.z},
        v1q = {q.x-p1.x, q.y-p1.y, q.z-p1.z},
     v12xv13 = {v12.y*v13.z-v12.z*v13.y,
                v12.z*v13.x-v12.x*v13.z,
                v12.x*v13.z-v12.z*v13.x};
   return (v1p.x*v12xv13.x + v1p.y*v12xv13.y + v1p.z*v12xv13.z) *
          (v1q.x*v12xv13.x + v1q.y*v12xv13.y + v1q.z*v12xv13.z) > 0;
}

// Point dans le parallépipède retenu
int dedans(POINT p)
{return p.x >= xmin && p.x <= xmax &&
        p.y >= ymin && p.y <= ymax && p.z >= zmin && p.z <= zmax;}


double distance(POINT p, POINT q){//TODO remove sqr
  return sqr((p.x-q.x)*(p.x-q.x)+
             (p.y-q.y)*(p.y-q.y)+
             (p.z-q.z)*(p.z-q.z));
}

double distanceMax(POINT a, POINT b, POINT c){
  double res1 = distance(a,b);
  double res2 = distance(a,c);
  double res3 = distance(c,b);
  double res = (res1<res2)?res2:res1;
  return (res<res3)?res3:res;
}

typedef struct {int n; POINT origine; POINT *point;
                unsigned *id;int32_t* phi;int32_t* theta;
                int nb_retenus;
                int err;} LECTURE_FICHIER;

LECTURE_FICHIER lecture_fichier(vector<POINT> vect, POINT origine){

  LECTURE_FICHIER res;//result
  res.err = TRUE;

  clock_t cpu = clock();

  int i, n; // nombre total de points
  POINT *point;  //Coordonnées (x,y,z) de chaque point de mesure

  n = vect.size();

  point = (POINT *) malloc(n * sizeof(POINT));
  res.id = (unsigned *) malloc((n+3) * sizeof(unsigned));

  // attention, dans la version avec fichier en param, on ne prend que les points dans une certaine zone
  for (i = 0; i < n; i++)
  {
    point[i] = vect[i];
  }

  /***** Sélection des points à retenir, calcul des azimuts et élévations ***/

  int nb_retenus = 0;
  for (i = 0; i < n; i++){
      res.id[nb_retenus++] = i;
  }
  res.id[nb_retenus] = nb_retenus;
  res.id[nb_retenus+1] = nb_retenus+1;
  res.id[nb_retenus+2] = nb_retenus+2;

  res.phi = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  res.theta = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  for (i = 0; i < nb_retenus; i++)
    calcule_phi_theta(origine, point[res.id[i]], res.phi+i, res.theta+i);

  printf("lecture d'un vector lu\n");
  printf("%d points lus en %f secondes\n", n, seconds(cpu));
  printf("%d nb_retenus\n",nb_retenus);

  res.n = n;
  res.origine = origine;
  res.point = point;
  res.nb_retenus = nb_retenus;
  res.err = FALSE;

  return res;
}

LECTURE_FICHIER lecture_fichier(const char file_path[]){
  int nb_lu;        //nb valeurs lues par fread
  FILE *fichier;
  unsigned taille_suppl; //information supplémentaire pour chaque point

  LECTURE_FICHIER res;//result
  res.err = TRUE;

  clock_t cpu = clock();

  fichier = fopen(file_path, "rb");

  int i, n;      // nombre total de points
  POINT origine; //Coordonnées de l'appareil de mesure
  POINT *point;  //Coordonnées (x,y,z) de chaque point de mesure

  nb_lu = fread(&origine.x, sizeof origine.x, 1, fichier);
  nb_lu = fread(&origine.y, sizeof origine.y, 1, fichier);
  nb_lu = fread(&origine.z, sizeof origine.z, 1, fichier);
  nb_lu = fread(&n, sizeof n, 1, fichier);
  nb_lu = fread(&taille_suppl, sizeof taille_suppl, 1, fichier);

  cout << "----------------" << endl << taille_suppl << endl << "--------------------" << endl;

  point = (POINT *) malloc(n * sizeof(POINT));
  res.id = (unsigned *) malloc((n+3) * sizeof(unsigned));

  int nb_retenus = 0;
  POINT temp;
  for (i = 0; i < n; i++)
  {
    nb_lu = fread(&temp.x, sizeof(temp.x), 1, fichier);
    nb_lu = fread(&temp.y, sizeof(temp.y), 1, fichier);
    nb_lu = fread(&temp.z, sizeof(temp.z), 1, fichier);

    //TODO pas la meilleur façon de faire, tableau ID surement inutile. Mieux comprendre l'algo pour supprimer le tableau d'ID
//    if (dedans(temp)){
      point[nb_retenus++] = temp;
//    }
    fseek(fichier, taille_suppl, SEEK_CUR);
  }
  if (!nb_lu) {return res;}
  fclose(fichier);

  /***** Sélection des points à retenir, calcul des azimuts et élévations ***/

  for (i = 0; i < nb_retenus; i++){
      res.id[i] = i;
  }
  res.id[nb_retenus]   = nb_retenus;
  res.id[nb_retenus+1] = nb_retenus+1;
  res.id[nb_retenus+2] = nb_retenus+2;

  res.phi = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  res.theta = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  for (i = 0; i < nb_retenus; i++)
    calcule_phi_theta(origine, point[res.id[i]], res.phi+i, res.theta+i);

  printf("%s lu\n", file_path);
  printf("%d points lus en %f secondes\n", n, seconds(cpu));
  printf("%d nb_retenus\n",nb_retenus);

  res.n = n;
  res.origine = origine;
  res.point = point;
  res.nb_retenus = nb_retenus;
  res.err = FALSE;

  return res;
}

typedef struct {LECTURE_FICHIER inputs;
                t_triangle * arbre;
                int nb_triangles_finaux;
                int nb_triangles;
                int err;} OUTPUT_TRIANGULATION;

OUTPUT_TRIANGULATION triangulation(LECTURE_FICHIER in, const double LIMITE){
  OUTPUT_TRIANGULATION res;
  res.inputs = in;
  res.err = TRUE;


  // pour que le code en dessous soit moins chargé
  int nb_retenus = in.nb_retenus;
  int32_t *phi = in.phi;
  int32_t *theta = in.theta;
  POINT* point = in.point;
  t_triangle* arbre = res.arbre;

  arbre = (t_triangle*) malloc(nb_triangles_par_point*nb_retenus*sizeof(t_triangle));

  if (!arbre)
  { printf("pas assez de memoire\n"); return res;}

  clock_t cpu = clock();
  phi[nb_retenus]   = INT_MIN/2; theta[nb_retenus]   = INT_MIN/2;
  phi[nb_retenus+1] = INT_MAX/2; theta[nb_retenus+1] = INT_MIN/2;
  phi[nb_retenus+2] = 0;         theta[nb_retenus+2] = INT_MAX/2;
  res.nb_triangles = construit_delaunay(nb_retenus, phi, theta, arbre);
  printf("arbre de %d triangles construits en %f secondes\n",
         res.nb_triangles, seconds(cpu));

  cpu = clock();
  int nb_triangles_finaux = 0;

  for (int i = 1; i < res.nb_triangles; i++)
  {
    if (arbre[i].final){
      if (arbre[i].p1 < nb_retenus && arbre[i].p2 < nb_retenus &&
          arbre[i].p3 < nb_retenus && distanceMax(point[arbre[i].p1],point[arbre[i].p2],point[arbre[i].p3])<LIMITE/**/){
        nb_triangles_finaux++;
      }
      else {
        arbre[i].final = FALSE;
      }
    }
   }

  printf("triangles finaux %d\n", nb_triangles_finaux);

  res.arbre = arbre;
  res.nb_triangles_finaux = nb_triangles_finaux;
  res.err = FALSE;
  return res;
}

void ecrire_fichier(const string output_path, OUTPUT_TRIANGULATION in, COULEUR couleurs[], int couleurID){
  clock_t cpu = clock();
  FILE *fichier = fopen(output_path.c_str(), "wb");
  fprintf(fichier, "ply\nformat binary_little_endian 1.0\nelement vertex %d\n",
          in.inputs.nb_retenus);
  fprintf(fichier, "property double x\nproperty double y\nproperty double z\n");
  fprintf(fichier, "element face %d\n", in.nb_triangles_finaux);
  fprintf(fichier, "property list uchar int vertex_index\n");
  fprintf(fichier, "property uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n");

  if(couleurID>11) {couleurID=11;}

  for (int i = 0; i < in.inputs.nb_retenus; i++){
      fwrite(&in.inputs.point[i].x,1,sizeof(double),fichier);
      fwrite(&in.inputs.point[i].y,1,sizeof(double),fichier);
      fwrite(&in.inputs.point[i].z,1,sizeof(double),fichier);
  }

  for (int i = 1; i < in.nb_triangles; i++)
    if (in.arbre[i].final)
      if (in.arbre[i].p1 >= in.inputs.nb_retenus || in.arbre[i].p2 >= in.inputs.nb_retenus || in.arbre[i].p3 >= in.inputs.nb_retenus)
        ; // triangle fictif
      else{
        fwrite(&nb_cote_triangle, 1, sizeof(unsigned char), fichier);
        fwrite(&in.arbre[i].p3, 1, sizeof(int), fichier);
        fwrite(&in.arbre[i].p2, 1, sizeof(int), fichier);
        fwrite(&in.arbre[i].p1, 1, sizeof(int), fichier);

         /* attribut une couleur en fonction des stations qui peuvent la voir */ // TODO FAIRE DE MANIERE PROCEDURALE (masque de bit)
         fwrite(&couleurs[couleurID].x, 1, sizeof(unsigned char), fichier);
         fwrite(&couleurs[couleurID].y, 1, sizeof(unsigned char), fichier);
         fwrite(&couleurs[couleurID].z, 1, sizeof(unsigned char), fichier);
        }
    else
      ; // triangle extérieur à la scène
  fprintf(fichier, "\n");

  fclose(fichier);

  printf("Fichier %s de %d points et %d triangles ecrit en %f secondes\n",
         output_path.c_str(), in.inputs.nb_retenus, in.nb_triangles_finaux, seconds(cpu));
}

void set_station_on_points(const OUTPUT_TRIANGULATION in,const POINT stations[], const int nb_stations){//TODO changer station en vecteur pour eviter d'avoir le nb_stations
  for (int i = 1; i < in.nb_triangles; i++){
    if (in.arbre[i].final){
      if (in.arbre[i].p1 >= in.inputs.nb_retenus || in.arbre[i].p2 >= in.inputs.nb_retenus || in.arbre[i].p3 >= in.inputs.nb_retenus)
        ; // triangle fictif
      else {
        for(int j=0;j<nb_stations;j++){ //TODO amelioration, on peut ignorer la station origine et directement mettre tous les points avec le masque sans passé par visible
          if(visible(in.inputs.origine,stations[j],
                      in.inputs.point[in.arbre[i].p3],
                      in.inputs.point[in.arbre[i].p2],
                      in.inputs.point[in.arbre[i].p1]))
          {
            in.inputs.point[in.arbre[i].p1].stations |= pow2(j);
            in.inputs.point[in.arbre[i].p2].stations |= pow2(j);
            in.inputs.point[in.arbre[i].p3].stations |= pow2(j);
          }
        }
      }
    }
  }
}

void get_points_by_allmask(map<unsigned int,vector<POINT>> &resultat, const OUTPUT_TRIANGULATION in){//in pour input
  for (int i = 0; i < in.inputs.nb_retenus; i++){
      resultat[in.inputs.point[i].stations].push_back(in.inputs.point[i]);
  }
}

void free_triangulation_partial(OUTPUT_TRIANGULATION out){
  free(out.inputs.id);
  free(out.inputs.phi);
  free(out.inputs.theta);
}

void free_triangulation_complete(OUTPUT_TRIANGULATION out){
  free_triangulation_partial(out);
  free(out.arbre);
  free(out.inputs.point);
}

template <typename T>
void clean_vector(vector<T> _v){
  _v.clear();
  _v.shrink_to_fit();
}

string binary(unsigned x){
    string s;
    do{
        s.push_back('0' + (x & 1));
    } while (x >>= 1);

    while(s.size()<15){
      s.push_back('0');
    }
    reverse(s.begin(), s.end());
    return s;
}

string addSlashEnd(string str){
  if(str.back() != '/')
    str += "/";
  return str;
}

vector<string> listFile(string path){
  path = addSlashEnd(path);

  vector<string> result;
  DIR *pDIR;
  struct dirent *entry;
  if( (pDIR=opendir(path.c_str())) ){
    while((entry = readdir(pDIR))){
      if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){
        result.push_back(path + entry->d_name);
      }
    }
    closedir(pDIR);
    return result;
  }
}

/*debug*/
void deb_sizes(map<unsigned int,vector<POINT>> _map){
  printf("\n nombre de vector : %d\n",_map.size());
  for(map<unsigned int,vector<POINT>>::iterator iter = _map.begin(); iter != _map.end(); ++iter){
    unsigned int k =  iter->first;

    printf("vect %d -> size : %d\n",k,_map[k].size());
  }
  printf("------------");
}

long long memory(){
  struct sysinfo memInfo;
  sysinfo (&memInfo);
  return memInfo.totalram-memInfo.freeram;
}

long long keepMaxMemory(long long old){
  long long mem = memory();
  if (old > mem)
    return old;
  return mem;
}

int main(int argc, char* argv[])
{
  auto memoryAtStart = memory();
  long long maxMemory = 0;
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();

  COULEUR couleurs[12];
  couleurs[0].x=0  ; couleurs[0].y=0  ; couleurs[0].z=255;
  couleurs[1].x=0  ; couleurs[1].y=255; couleurs[1].z=0  ;
  couleurs[2].x=255; couleurs[2].y=0  ; couleurs[2].z=0  ;
  couleurs[3].x=0  ; couleurs[3].y=255; couleurs[3].z=255;
  couleurs[4].x=255; couleurs[4].y=0  ; couleurs[4].z=255;
  couleurs[5].x=255; couleurs[5].y=255; couleurs[5].z=0  ;
  couleurs[6].x=0  ; couleurs[0].y=0  ; couleurs[0].z=128;
  couleurs[7].x=0  ; couleurs[1].y=128; couleurs[1].z=0  ;
  couleurs[8].x=128; couleurs[2].y=0  ; couleurs[2].z=0  ;
  couleurs[9].x=0  ; couleurs[3].y=128; couleurs[3].z=128;
  couleurs[10].x=128; couleurs[4].y=0  ; couleurs[4].z=128;
  couleurs[11].x=128; couleurs[5].y=128; couleurs[5].z=0  ;

  if (!strcmp(argv[1],"-h"))
  {
    printf("usage: %s precisionNonMerged precisionMerged inputFolder outputFolder [-o]\nIf folder already exit, do nothing. If folder doesn't exist, it create it", argv[0]);
    return EXIT_FAILURE;
  }

  bool saveOriginal = (argc == 6 && !strcmp(argv[5],"-o"));

  /* récupérer les positions des stations */
  string inputPath = argv[3];
  vector<string> inputs = listFile(inputPath);

  const int nb_stations = inputs.size();
  POINT stations[nb_stations];//TODO opti struct station qui contient la position, et son file_path
  for(int i=0;i<nb_stations;i++){
    stations[i] = extrait_station(inputs[i].c_str());
  }

  double LIMITE1 = atof(argv[1]);
  double LIMITE2 = atof(argv[2]);
  string output_name(argv[4]);
  output_name = addSlashEnd(output_name);
  string output_original = output_name+"original/";
  map<unsigned int,vector<POINT>> points;

  /******************* lecture du fichier de données ************************/
  /***** Sélection des points à retenir, calcul des azimuts et élévations ***/
  for(int i=0;i<nb_stations;i++){
    //dans le cas où l'utilisateur veut les triangulation original, tout est fait ici
    if(saveOriginal){
      LECTURE_FICHIER input_temp = lecture_fichier(inputs[i].c_str());
      OUTPUT_TRIANGULATION out_temp = triangulation(input_temp, LIMITE1);
      ecrire_fichier(output_original+to_string(i)+".ply",out_temp,couleurs,i);
      free_triangulation_partial(out_temp);
      free(out_temp.arbre);
    }

    LECTURE_FICHIER input_tri = lecture_fichier(inputs[i].c_str());
    OUTPUT_TRIANGULATION out_tri = triangulation(input_tri, 100000000);// 100000000 pour ignorer cette règle lors des première triangulation

    if (input_tri.err == TRUE) {return EXIT_FAILURE;}
    if (out_tri.err   == TRUE) {return EXIT_FAILURE;}

    // pour les metriques
    maxMemory = keepMaxMemory(maxMemory);

    free_triangulation_partial(out_tri);
    set_station_on_points(out_tri,stations,nb_stations);
    free(out_tri.arbre);
    get_points_by_allmask(points,out_tri);
  }

  /* TRIANGULATION GENERALISÉ */
  for(map<unsigned int,vector<POINT>>::iterator iter = points.begin(); iter != points.end(); ++iter){
    unsigned int k =  iter->first;
    if(k == 0){ //TODO pourquoi il y a (un)des points qui ne sont vu par aucune station ?
      continue;
    }

    // récupère le premier bit à 1 pour définir la station utilisé
    // dans le cas d'un algo  mergé, retourne le premier bit à 1
    // (vu qu'il faut triangulé en fonction d'une des stations)
    int station = -1;
    while(!(k&pow2(++station))) {}

    // debug
//    for(int i=0;i<points.size();i++){
//      if (points[i].size()==0) {continue;}
//      cout << points[i].size() << endl;
//    }

    // ignore les triangulation qui ont trop peu de points
    if(points[k].size()<100){ // le 100 est arbitraire, on pourrait surement mettre un chiffre plus bas.
      cout << output_name+binary(k)+".ply est ignoré car contient "<< points[k].size() <<"<100 points" << endl;
      continue;
    }

    // c'est une station unique si seulement 1 bit existe
    // sinon c'est une des stations selectionné
    LECTURE_FICHIER lf_station1 = lecture_fichier(points[k],stations[station]);
    if (lf_station1.err == TRUE){return EXIT_FAILURE;}
    // si une seule station LIMITE 1
    // si plusieur station, une des stations donc LIMITE2
    OUTPUT_TRIANGULATION ot_station1 = triangulation(lf_station1,(isPow2(k))?LIMITE1:LIMITE2);
    maxMemory = keepMaxMemory(maxMemory);
    if (ot_station1.err == TRUE){return EXIT_FAILURE;}
    ecrire_fichier(output_name+binary(k)+".ply",ot_station1,couleurs,k-1);//TODO afficher le nombre en binaire
    clean_vector(points[k]);
    free_triangulation_complete(ot_station1);
  }

  //affichage de metriques
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  cout << "Temps total de l'execution : " << elapsed_seconds.count() << "s" <<endl;
  cout << "Usage maximum de la maymoire : " << (maxMemory-memoryAtStart)/(1000000.0) << " MB" << endl;

  return EXIT_SUCCESS;
}

