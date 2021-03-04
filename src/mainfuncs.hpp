#ifndef MAINFUNCS
#define MAINFUNCS

#include "main.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"


/* Classes and structures */
struct overrides_ {
  //io settings

  bool input = false;
  bool charges_input = false;
  bool types_input = false;
  bool output = false;
  bool md_output = false;

  //spatial settings

  bool layers = false;
  bool diameter = false;
  bool distance = false;

  //potential settings

  bool ljsteps = false;
  bool ljstepsize = false;
  bool cdsteps = false;
  bool cdstepsize = false;
  bool ljmin = false;
  bool cdmin = false;

  bool cdcut = false;
  bool ljcut = false;

  // MD settings
  bool MD_walltime = false;
  bool MD_cores = false;
  bool MD_numMolperDim = false;
  bool MD_dielectric = false;
  bool MD_LJepsilon = false;
  bool MD_LJcutoff = false;
  bool MD_CDcutoff = false;
  bool MD_timestep = false;
  bool MD_runtime = false;
};


/* Functions */
// std::string hash_output(collagenFibril &fib);
int process_arg(std::string strarg, int* errcodes, bool dashes, collagenFibril &fib);
int read_config_file(std::string path, collagenFibril &fib);
int read_args(int argc, char const *argv[], collagenFibril &fib);
int parse_all_args(int argc, char const *argv[], collagenFibril &fib);
void readAtomInfos(collagenFibril &fib);
void programInfo();
void printOptions(collagenFibril fib);
void print_help();

#endif
