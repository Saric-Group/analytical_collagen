#ifndef MAINFUNCS
#define MAINFUNCS

#include "main.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"


/* Classes and structures */
struct overrides_ {
  //io settings

  bool input_override = false;
  bool output_override = false;

  //spatial settings

  bool layers_override = false;
  bool diameter_override = false;
  bool distance_override = false;

  //potential settings

  bool ljsteps_override = false;
  bool ljstepsize_override = false;
  bool cdsteps_override = false;
  bool cdstepsize_override = false;
  bool ljmin_override = false;
  bool cdmin_override = false;

  bool cdcut_override = false;
  bool ljcut_override = false;
};


/* Functions */
std::string hash_output(collagenFibril &fib);
int process_arg(std::string strarg, int* errcodes, bool dashes, collagenFibril &fib);
int read_config_file(std::string path, collagenFibril &fib);
int read_args(int argc, char const *argv[], collagenFibril &fib);
int parse_all_args(int argc, char const *argv[], collagenFibril &fib);
void print_help();

#endif
