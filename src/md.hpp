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
  double x;
  double y;
  double z;
  double length = 0.;

  /* Constructors */
  position(double x_ = 0., double y_ = 0., double z_ = 0.)
  {
    x = x_;
    y = y_;
    z = z_;
    calcLength();
  }

  /* Functions */
  void calcLength();
  void normalize();
  void printPos();
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
    for (int i = 1; i <= numPoints; i++) {
      x = i * shift - length / 2.0;
      for (int j = 1; j <= numPoints; j++) {
        y = j * shift - length / 2.0;
        for (int k = 1; k <= numPoints; k++) {
          z = k * shift - length / 2.0;
          position pos{x, y, z};
          gridPoints.push_back(pos);
        }
      }
    }
  }

  /* Functions */
};


/* Functions */
double dot(position a, position b);
position normal(position a, position b);
position getCapsuleA(position centerPoint, double length, double phi,
                          double theta);
position getCapsuleB(position centerPoint, double length, double phi,
                          double theta);
position closestPtOnLineSegment(position a, position b, position p);
bool capsuleOverlap(int a, int b, collagenFibril fib,
                  std::vector<position> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem);
bool checkOverlap(int i, int L, std::vector<position> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem,
                  collagenFibril fib);

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
