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
#include <sys/stat.h>       /* chmod */


typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;


/* Classes and structures */
struct filePaths_ {
  std::string mainpath = "";
  std::string inputpath = "./molecule/charges_distribution_36";
  std::string outputpath = "./energy_min.dat";
  std::string file_extension = ".dat";
  std::string csv_extension = ".csv";
  std::string configpath = "";

  std::string md_outputpath = "./md/";
};

struct flags_ {
  bool consoleOutput = false;
  bool set_output = false; // @Joel: what is this used for?
  bool csv_output = false;
  bool charge_hashed_outputs = false;
  bool xyz_outputs = false;
  bool printAtomInfo = false;
  bool printMoleculeInfo = false;
  bool originalOutput = false;
  bool mdOutput = false;
  bool measureTime = false;
  bool development = false;
};

#endif
