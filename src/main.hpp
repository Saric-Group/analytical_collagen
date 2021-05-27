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
#include <stdio.h>


typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;


/* Classes and structures */
struct filePaths_ {
  std::string mainpath;
  std::string inputpath;
  std::string outputpath;
  std::string file_extension = ".dat";
  std::string csv_extension = ".csv";
  std::string configpath;

  std::string md_outputpath;
};

struct flags_ {
  bool help = false;
  bool consoleOutput = false;
  bool debugInfo = false;
  bool printAtomInfo = false;
  bool printMoleculeInfo = false;
  bool layermodel = false;
  bool minimize = false;
  bool lammps = false;
  bool md = false;
  bool measureTime = false;
  bool development = false;
};

struct dev_ {
  bool random;
  int seed;
  double factor;
  int samples;
  int pos, neg;
};

#endif
