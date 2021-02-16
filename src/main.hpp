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


/* Classes and structures */
struct filePaths_ {
  std::string inputpath = "./charge_distribution";
  std::string outputpath = "./energy_min.dat";
  std::string file_extension = ".dat";
  std::string configpath = "";
};

struct flags_ {
  bool set_output = false;
  bool charge_hashed_outputs = true;
  bool xyz_outputs = false;
};

#endif
