#ifndef COLLMOL
#define COLLMOL

#include "main.hpp"


/* Classes and structures */
struct collagenMolecule {
  int numAtoms = 0;
  int numTypes = 0;
  int numPos = 0;
  int numNeg = 0;

  double distanceAtoms = 0.255;
  double diameterAtom = 1.12;
  double chargeAtom = 22.4;
  double length;
  double totalCharge = 0.;

  std::vector<double> charges, nonZeroChargesIndex;
  std::vector<int> atomTypes;

  const int maxNumTypes = 20;
  int typeCharge [20] = {0,   // 01 Glutamine GLN
                         0,   // 02 Methionine MET
                         0,   // 03 Serine SER
                         0,   // 04 Tyrosine TYR
                         0,   // 05 Glycine GLY
                         -1,   // 06 Aspartic Acid ASP
                         -1,   // 07 Glutamic Acid GLU
                         1,   // 08 Lysine LYS
                         0,   // 09 Alanine ALA
                         0,   // 10 Valine VAL
                         0,   // 11 Proline PRO
                         1,   // 12 Arginine ARG
                         0,   // 13 Leucine LEU
                         0,   // 14 (Hydroxyproline) HYP
                         0,   // 15 Phenylalanine PHE
                         0,   // 16 Asparagine ASN
                         0,   // 17 Threonine THR
                         0,   // 18 ? LYZ
                         0,   // 19 Histidine HIS
                         0   // 20 Isoleucine ILE
                       };

  /* Constructors */

  /* Functions */
  void typesFromCharges();
  void countCharges();
  void chargesFromTypes();
  void readCharges(std::string &file);
  void readCharges(std::vector<double> vec);
  void readTypes(std::string &file);
  void readTypes(std::vector<int> vec);

  void printAtoms();
  void printMoleculeInfo();
};

#endif
