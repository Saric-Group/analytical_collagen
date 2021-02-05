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
struct parameters_ {
  double lj_min = 0.01;
  double lj_stepsize = 0.01;
  int lj_steps = 50;

  double cd_min = 10.0;
  double cd_stepsize = 10.0;
  int cd_steps = 20;

  double cd_cutoff = 2.0;     /* cutoff of cd potential */
  double lj_cutoff = 2.0;     /* cutoff of lj potential */
  double max_cutoff = std::max(cd_cutoff, lj_cutoff);
};

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
