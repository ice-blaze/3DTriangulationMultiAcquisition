# 3DTriangulationMultiAcquisition



##### Don't forget (fridge)
/usr/include/CGAL/Interval_nt.h


##### TODO
* extraire les points des deux stations
* structure de ce qu'on a fait, pour le rapport
* triangulation cgal
* 2 triangulation séparé dans un premier coup et ensuite la triangulation a merger

##### Comparatif PERFORMANCE sur fichier à 6 millions de points
* CGAL triangulation normale
  * 14m11s
* CGAL triangulation Delaunay
  * 13m55s
* CGAL KD
  * 2m24s
* KD Perso
  * 8s
* Triangulation Perso
  * -

##### Structure du Rapport
* Introduction
* Description du problème
* Développement de l'algorithme
  * Fonctionnement du principe
  * Structures de données utilisées
  * Technologie(c/c++/CGAL)
  * Perfomances (complexité, speedups)
* Conclusion (synthèse des resultats futur)
