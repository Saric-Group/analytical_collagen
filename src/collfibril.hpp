#ifndef COLLFIBRIL
#define COLLFIBRIL

#include "main.hpp"
#include "parse.hpp"
#include "funcs.hpp"
#include "collmol.hpp"


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

  /* Functions */
  void calcCutoff() {
    max_cutoff = std::max(cd_cutoff, lj_cutoff);
  }
};

struct parametersMD_ {
  bool outputScript = false;
  bool outputTopology = false;
  bool outputLAMMPSinput = false;
  bool scriptbuild = false;
  bool rigid = false;

  int walltime = 24;
  int cores = 1;

  int numMolperDim = 5;

  double kAngle = 50.0;
  double kAngle_start = 50.0;
  double kAngle_inc = 50.0;
  double kAngle_end = 50.0;

  double dielectric = 10.0;
  double dielectric_start = 10.0;
  double dielectric_inc = 10.0;
  double dielectric_end = 100.0;

  double LJepsilon = 0.01;
  double LJepsilon_start = 0.01;
  double LJepsilon_inc = 0.01;
  double LJepsilon_end = 0.5;

  double cd_cutoff = 5.0;
  double lj_cutoff = 5.0;

  double timestep = 0.002;
  int runtime = 6000001;

  std::string lmp_mpi = "~/Scratch/lammps-29Oct20/src/lmp_mpi";
};

struct collagenFibril {
  collagenMolecule mol;
  parameters_ parameters;
  parametersMD_ parametersMD;
  int layers = 2;
  double latGap = mol.diameterAtom;
  double radGap = mol.diameterAtom;
  double offset = mol.diameterAtom;
  double energy = 1e16;

  /* Constructors */

  /* Old functions */
  double distance(double pos, double first, int n, double lat_gap);

  double factorLJ(double d);
  double LJ_per_mol(double pos, double dx, double ref, double lat_gap);
  double LJ_layer(double pos, double dx, int layer, double offset,
                  double lat_gap, double box);
  double compLJ(double lat_gap, double rad_gap, double offset);

  double factorCD(double q1, double q2, double d);
  double CD_per_mol(double pos, double q1, double dx, double ref, double lat_gap);
  double CD_layer(double pos, double q1, double dx, int layer, double offset,
                  double lat_gap, double box);
  double compCD(double lat_gap, double rad_gap, double offset);

  std::string hash_output();
  int produce_xyz(int numatoms, int rows, double distatom, double latgap, double radgap, double offset, int id[], int type[], double charges[], double xpos[], double ypos[], double zpos[]);
  int output_xyz(std::string filename, int rows, double latgap, double radgap, double offset);
  void header(std::string &file);
  void singleEmin();

  /* New functions */
  void minimizeEnergy();
  void writeXYZ();
  void printMoleculeInfo();
};

#endif
