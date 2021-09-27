#ifndef MD
#define MD

#include "main.hpp"
#include "layermodel.hpp"
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

struct vec3d {
  double x, y, z;
  double length = 0.;

  /* Constructors */
  vec3d(double x_ = 0., double y_ = 0., double z_ = 0.)
  {
    x = x_;
    y = y_;
    z = z_;
    calcLength();
  }

  /* Functions */
  void calcLength();
  void normalize();

  /* Operators */
  vec3d operator+(const vec3d &a) const
  {
    return vec3d(x + a.x, y + a.y, z + a.z);
  }
  vec3d operator-(const vec3d &a) const
  {
    return vec3d(x - a.x, y - a.y, z - a.z);
  }
  vec3d operator*(const double &k) const
  {
    return vec3d(k * x, k * y, k * z);
  }
  double operator*(const vec3d &a) const
  {
    return a.x * x + a.y * y + a.z * z;
  }
  vec3d operator*=(const double &k)
  {
    x *= k;
    y *= k;
    z *= k;
    return *this;
  }
  vec3d operator/(const double &k)
  {
    return vec3d(x / k, y / k, z / k);
  }
  vec3d operator/=(const double &k)
  {
    x /= k;
    y /= k;
    z /= k;
    return *this;
  }
};
vec3d operator*(const double &k, const vec3d &a);
std::ostream& operator<<(std::ostream &os, const vec3d &a);


struct cubeGrid {
  double length;    // length of cube side
  double shift;     // distance of nearest neighbours
  int numPoints;    // points per side, total #points = numPoints^3
  std::vector<vec3d> gridPoints;


  /* Constructors */
  cubeGrid(double length_ = 100, int numPoints_ = 1)
  {
    length = length_;
    numPoints = numPoints_;
    shift = length / (numPoints + 1);
    double x, y, z;
    for (int i = 1; i <= numPoints; i++) {
      x = i * shift - length / 2.0;
      for (int j = 1; j <= numPoints; j++) {
        y = j * shift - length / 2.0;
        for (int k = 1; k <= numPoints; k++) {
          z = k * shift - length / 2.0;
          vec3d pos(x, y, z);
          gridPoints.push_back(pos);
        }
      }
    }
  }

  /* Functions */
};


/* Functions */
vec3d crossproduct(vec3d a, vec3d b);
vec3d getLinePtA(vec3d centerPoint, double length, double phi, double theta);
vec3d getLinePtB(vec3d centerPoint, double length, double phi, double theta);
double shortestDistBetweenLines(vec3d a1, vec3d a2, vec3d b1, vec3d b2);
bool moleculesOverlap(int a, int b, collagenMolecule mol,
                  std::vector<vec3d> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem);
bool checkOverlap(int i, int L, std::vector<vec3d> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem,
                  collagenMolecule mol);

void genTopologyZero(collagenMolecule mol, int L);
void printScriptVars(FILE *outf, md_var var, int tabs);
std::vector<md_var> collectMDvars(collagenMolecule mol);
std::string getDir(std::vector<md_var> md_vars, bool pvalues = false);
std::string getCargs(std::vector<md_var> md_vars);
void genInSim(collagenMolecule mol);
void genQsub(collagenMolecule mol);
void genBashScript(collagenMolecule mol);
void createLAMMPSfiles(collagenMolecule mol);
#endif
