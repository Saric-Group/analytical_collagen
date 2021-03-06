#ifndef MD
#define MD

#include "main.hpp"
#include "collfibril.hpp"
#include "random.hpp"


/* Classes and structures */
struct md_var {
  std::string abr;
  int size;
  int prec = 3;
  std::vector<void *> vars;
  void *var;

  /* Constructors */
  md_var(std::string abr_, int size_, double *var_, std::vector<double *> vars_, int prec_) {
    abr = abr_;
    size = size_;
    prec = prec_;
    var = var_;
    for (int i = 0; i < (int) vars_.size(); i++) {
      vars.push_back(vars_[i]);
    }
  }
  md_var(std::string abr_, int size_, int *var_, std::vector<int *> vars_) {
    abr = abr_;
    size = size_;
    var = var_;
    for (int i = 0; i < (int) vars_.size(); i++) {
      vars.push_back(vars_[i]);
    }
  }
};

struct position {
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

struct cubeGrid {
  double length;    // length of cube side
  int numPoints;    // points per side, total #points = numPoints^3
  std::vector<position> gridPoints;


  /* Constructors */
  cubeGrid(double length_, int numPoints_)
  {
    length = length_;
    numPoints = numPoints_;
    double shift = length / (numPoints + 1);
    double x, y, z;
    position pos;
    for (int i = 1; i <= numPoints; i++) {
      x = i * shift - length / 2.0;
      for (int j = 1; j <= numPoints; j++) {
        y = j * shift - length / 2.0;
        for (int k = 1; k <= numPoints; k++) {
          z = k * shift - length / 2.0;
          pos.x = x;
          pos.y = y;
          pos.z = z;
          gridPoints.push_back(pos);
        }
      }
    }
  }

  /* Functions */
};


/* Functions */
void genTopologyZero(collagenFibril fib, int L);
void printScriptVars(FILE *outf, md_var var, int tabs);
std::vector<md_var> collectMDvars(collagenFibril fib);
std::string getDir(std::vector<md_var> md_vars, bool pvalues = false);
std::string getCargs(std::vector<md_var> md_vars);
void genInSim(collagenFibril fib);
void genQsub(collagenFibril fib);
void genBashScript(collagenFibril fib);
void createLAMMPSfiles(collagenFibril fib);
#endif
