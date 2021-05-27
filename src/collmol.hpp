#ifndef COLLMOL
#define COLLMOL

#include "main.hpp"


/* Classes and structures */
struct parametersMD_ {
  bool script = false;
  bool topology = false;
  std::string topFile = "topology.0time";
  bool input = false;
  std::string inputFile = "in.sim";
  std::string logFile = "log.sim";
  std::string dumpFile = "out_sim.xyz";
  bool scriptbuild = false;
  bool rigid = false;

  int walltime = 24;
  int cores = 1;

  int numMolperDim = 5;
  bool random = false;
  double phi = 0.0;
  double theta = 0.3 * M_PI;

  double kAngle = 50.0;
  double kAngle_start = 50.0;
  double kAngle_inc = 50.0;
  double kAngle_end = 50.0;

  double dielectric = 10.0;
  double dielectric_start = 10.0;
  double dielectric_inc = 10.0;
  double dielectric_end = 100.0;
  double cd_cutoff = 5.0;

  double LJepsilon = 0.01;
  double LJepsilon_start = 0.01;
  double LJepsilon_inc = 0.01;
  double LJepsilon_end = 0.5;
  double lj_cutoff = 5.0;

  double timestep = 0.002;
  int runtime = 6000001;

  std::string lmp_mpi = "~/Scratch/lammps-29Oct20/src/lmp_mpi";
};

struct collagenMolecule {
  parametersMD_ parametersMD;

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

  int maxNumTypes = 20;
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
  void moleculeToFile(std::string &file);
  void chargesToFile(std::string &file, int k = 0);
  void addOvitoHeaderToChargeFile(std::string &file);

  std::vector<double> smoothen(int delta);
  std::vector<double> binNormalize(std::vector<double> vec, int n);

};

#endif
