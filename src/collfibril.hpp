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
};

#endif
