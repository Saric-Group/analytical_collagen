#ifndef MAIN
#define MAIN

#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */
#include <vector>           /* vector, ... */
#include <sstream>
#include <map>
#include <algorithm>        /* std::remove */
#include <iomanip>          /* std::setw() */


/* Classes and structures */
struct filePaths_ {
  std::string inputpath = "./charge_distribution";
  std::string charges_inputpath = "./molecule/charges_distribution_36";
  std::string types_inputpath = "./molecule/atom_types_1054";
  std::string outputpath = "./energy_min.dat";
  std::string file_extension = ".dat";
  std::string configpath = "";

  std::string md_outputpath = "./md/";
};

struct flags_ {
  bool set_output = false; // @Joel: what is this used for?
  bool charge_hashed_outputs = false;
  bool xyz_outputs = false;
  bool readCharges = false;
  bool readTypes = false;
  bool printAtomInfo = false;
  bool printMoleculeInfo = false;
  bool annesOutput = false;
  bool measureTime = false;
};

#endif
