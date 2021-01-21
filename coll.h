#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */
#include <vector>           /* vector, ... */

using namespace std;

/* Functions */
void readAtoms(string &file);
void readTypes(string &file);

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

void header(string &file);
void singleEmin(string &file, int layers);
void multipleEmin(string &file, int layers, int numbers);
/* map functions */
double energy(double lat_gap, double rad_gap, double offset, int layer1,
              int layer2, double lj, double cd, bool lj_on);
double intra_layer_energy(double lat_gap, double rad_gap, double offset,
                          double lj, double cd, bool lj_on);
int read_args(int, char *);
