#ifndef COLLMOL
#define COLLMOL

#include "main.hpp"


/* Classes and structures */
struct collagenMolecule {
  int numAtoms = 0;
  double distanceAtoms = 0.285;
  double diameterAtom = 1.12;
  double length;
  std::vector<double> charges, nonZeroChargesIndex;
  std::vector<int> atomType;

  /* Constructors */

  /* Functions */
  void readAtoms(std::string &file);
  void readAtoms(std::vector<double> vec);
};

#endif
