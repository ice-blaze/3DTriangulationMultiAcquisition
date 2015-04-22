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

// Pour passer de coordonnées de type double en int32_t
#define PRECISION 10000000

// Les points de mesure à prendre en considération seront à l'intérieur d'un
// parallélipipède compris entre [(xmin, ymin, zmin), (xmax, ymax, zmax)]
//double xmin = -1.0, xmax = 3.0, ymin = -2.0, ymax = 3.0, zmin = -1.0, zmax = 3.0;
double xmin = -10.0, xmax = 30.0,
       ymin = -20.0, ymax = 30.0,
       zmin = -10.0, zmax = 30.0;

typedef struct {double x; double y; double z;} POINT;

const int FALSE = 0;
const int TRUE = 1;

// L'arborescence a environ 9 triangles par point; 10 pour être confortable
const int nb_triangles_par_point = 10;
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
{ int* permutation;
  int i, nb_triangles;

  permutation = (int*) malloc( n*sizeof(int));
  for (i = 0; i < n; i++) permutation[i] = i;
  for (i = 0; i < n-1; i++)
    transpose(&permutation[i], &permutation[unif(i, n-1)]);

  arbre[0] =(t_triangle){n, n+1, n+2, -1, -1, -1, TRUE};
  nb_triangles = 1;

  for (i = 0; i < n; i++)
  {
     if (!insere_point(permutation[i], x, y, &nb_triangles, arbre))
       printf("%de point pas insere; %d\n", i, nb_triangles);
     if (nb_triangles >= nb_triangles_par_point*n)
     {
       printf("Arbre mal dimensionne; augmenter nb_triangles_par_point\n");
       return 0;
     }
     if (i%(n/100) == 0)
     { printf("\b\b\b%2d%%", 100*i/n); fflush(stdout);}
  }
  printf("\n");

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
POINT extrait_station(char chemin_ply[]){

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

int main(int argc, char* argv[])
{

  clock_t cpu;

  int nb_triangles; //nb triangles de la triangulation
  int nb_retenus;   //nb points dans le parallélipipède
  int nb_lu;        //nb valeurs lues par fread
  int nb_triangles_finaux; //nb triangles à afficher
  int32_t* phi;            //azimuth de chaque point
  int32_t* theta;          //élévation de chaque point
  unsigned *id;            //identificateur des points retenus
  t_triangle * arbre;      //noeuds de l'arborescence de la triangulation
  FILE *fichier;

  if (!strcmp(argv[1],"-h"))
  {
    printf("usage: %s infile1.bin [infile2.bin ... infileN.bin] outfile.ply\n", argv[0]);
    return EXIT_FAILURE;
  }

  cpu = clock();

  /* récupérer les positions des stations */
  POINT stations[argc-2];
  for(int i=1;i<argc-1;i++){
    stations[i] = extrait_station(argv[i]);
    printf("%g %g %g\n",stations[i].x,stations[i].y,stations[i].z);
  }

  /******************* lecture du fichier de données ************************/
  fichier = fopen(argv[1], "rb");

  unsigned i, n; // nombre total de points
  unsigned taille_suppl; //information supplémentaire pour chaque point
  unsigned char nb_cote_triangle = 3; //Pour l'écriture de 3 en binaire
  POINT origine; //Coordonnées de l'appareil de mesure
  POINT *point;  //Coordonnées (x,y,z) de chaque point de mesure

  nb_lu = fread(&origine, sizeof origine, 1, fichier);
  nb_lu = fread(&n, sizeof n, 1, fichier);
  nb_lu = fread(&taille_suppl, sizeof taille_suppl, 1, fichier);

  point = (POINT *) malloc(n * sizeof(POINT));
  id = (unsigned *) malloc((n+3) * sizeof(unsigned));

  for (i = 0; i < n; i++)
  {  nb_lu = fread(point + i, sizeof(POINT), 1, fichier);
     fseek(fichier, taille_suppl, SEEK_CUR);
  }
  if (!nb_lu) return EXIT_FAILURE;
  fclose(fichier);

  /***** Sélection des points à retenir, calcul des azimuts et élévations ***/
  nb_retenus = 0;
  for (i = 0; i < n; i++)
    if (dedans(point[i]))
      id[nb_retenus++] = i;
  id[nb_retenus] = nb_retenus;
  id[nb_retenus+1] = nb_retenus+1;
  id[nb_retenus+2] = nb_retenus+2;

  phi = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  theta = (int32_t *) malloc((nb_retenus+3) * sizeof(int32_t));
  for (i = 0; i < nb_retenus; i++)
    calcule_phi_theta(origine, point[id[i]], phi+i, theta+i);

  printf("%s lu\n", argv[1]);
  printf("%d points lus en %f secondes\n", n, seconds(cpu));
  printf("%d nb_retenus\n",nb_retenus);

  /*************** Construction de la triangulation  ***************/

  arbre = (t_triangle*) malloc(nb_triangles_par_point*nb_retenus*sizeof(t_triangle));

  if (!arbre)
  { printf("pas assez de memoire\n"); return EXIT_FAILURE;}

  cpu = clock();
  phi[nb_retenus]   = INT_MIN/2; theta[nb_retenus]   = INT_MIN/2;
  phi[nb_retenus+1] = INT_MAX/2; theta[nb_retenus+1] = INT_MIN/2;
  phi[nb_retenus+2] = 0;         theta[nb_retenus+2] = INT_MAX/2;
  nb_triangles = construit_delaunay(nb_retenus, phi, theta, arbre);
  printf("arbre de %d triangles construits en %f secondes\n",
         nb_triangles, seconds(cpu));

  cpu = clock();
  nb_triangles_finaux = 0;
  for (i = 1; i < nb_triangles; i++)
  {
    if (arbre[i].final)
      if (arbre[i].p1 < nb_retenus && arbre[i].p2 < nb_retenus &&
          arbre[i].p3 < nb_retenus)
        nb_triangles_finaux++;
   }

  printf("triangles finaux %d\n", nb_triangles_finaux);

  /***********Écriture de la triangulation au format ply ************/
  fichier = fopen(argv[argc-1], "wb");
  fprintf(fichier, "ply\nformat binary_little_endian 1.0\nelement vertex %d\n",
          nb_retenus);
  fprintf(fichier, "property double x\nproperty double y\nproperty double z\n");
  fprintf(fichier, "element face %d\n", nb_triangles_finaux);
  fprintf(fichier, "property list uchar int vertex_index\nend_header\n");

  for (i = 0; i < nb_retenus; i++)
    fwrite(point + id[i], 1, sizeof(POINT), fichier);

  for (i = 1; i < nb_triangles; i++)
    if (arbre[i].final)
      if (arbre[i].p1 >= n || arbre[i].p2 >= n || arbre[i].p3 >= n )
        ; // triangle fictif
      else
        {fwrite(&nb_cote_triangle, 1, sizeof(unsigned char), fichier);
         fwrite(&arbre[i].p3, 1, sizeof(int), fichier);
         fwrite(&arbre[i].p2, 1, sizeof(int), fichier);
         fwrite(&arbre[i].p1, 1, sizeof(int), fichier);
        }
    else
      ; // triangle extérieur à la scène
  fprintf(fichier, "\n");

  fclose(fichier);
  printf("Fichier %s de %d points et %d triangles ecrit en %f secondes\n",
         argv[2], nb_retenus, nb_triangles_finaux, seconds(cpu));
  free(arbre);
  return EXIT_SUCCESS;
 }

