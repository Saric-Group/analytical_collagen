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






/*****************************************************************************/

/* Reorganize below functions */
// void readAtoms(string &file);
void readTypes(std::string &file);

double distance(double pos, double first, int n, double lat_gap);

double factorLJ(double d);
double LJ_per_mol(double pos, double dx, double ref, double lat_gap);
double newLJ_per_mol(double pos, double dx, double ref, double lat_gap);
double LJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap);
double newLJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap);
double compLJ(double lat_gap, double rad_gap, double offset, int layers);
double newLJ(double lat_gap, double rad_gap, double offset, int layers);

double factorCD(double q1, double q2, double d);
double CD_per_mol(double pos, double q1, double dx, double ref,
                  double lat_gap);
double CD_layer(double pos, double q1, double dx, int layer, double offset,
                double lat_gap);
double compCD(double lat_gap, double rad_gap, double offset, int layers);

void header(std::string &file);
void singleEmin(std::string &file, int layers, collagenFibril fib);
void multipleEmin(std::string &file, int layers, int numbers);
/* map functions */
double energy(double lat_gap, double rad_gap, double offset, int layer1,
              int layer2, double lj, double cd, bool lj_on);
double intra_layer_energy(double lat_gap, double rad_gap, double offset,
                          double lj, double cd, bool lj_on);
int parse_all_args(int, char const **);
int read_args(int, char const **);
int process_arg(std::string, int*, bool);
void print_help();
int read_config_file(std::string);


#endif
