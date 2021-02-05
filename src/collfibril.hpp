#ifndef COLLFIBRIL
#define COLLFIBRIL

#include "main.hpp"
#include "collmol.hpp"


/* Classes and structures */
struct collagenFibril {
  collagenMolecule mol;
  int layers = 0;
  double latGap = mol.diameterAtom;
  double radGap = mol.diameterAtom;
  double offset = mol.diameterAtom;
  double energy = 1e16;

  /* Constructors */

  /* Old functions */
  double distance(double pos, double first, int n, double lat_gap);

  double factorLJ(double d);
  double LJ_per_mol(double pos, double dx, double ref, double lat_gap);
  double LJ_layer(double pos, double dx, int layer, double offset,
                  double lat_gap, double box);
  double compLJ(double lat_gap, double rad_gap, double offset);

  double factorCD(double q1, double q2, double d);
  double CD_per_mol(double pos, double q1, double dx, double ref, double lat_gap);
  double CD_layer(double pos, double q1, double dx, int layer, double offset,
                  double lat_gap, double box);
  double compCD(double lat_gap, double rad_gap, double offset);

  void singleEmin();

  /* New functions */
  void minimizeEnergy();
  void writeXYZ();
};

#endif
