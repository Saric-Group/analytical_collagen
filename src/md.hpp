#ifndef MD
#define MD

#include "main.hpp"
#include "collfibril.hpp"
#include "random.hpp"


/* Classes and structures */
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
void genInSim(collagenFibril fib);
void genQsub(collagenFibril fib);
void genBashScript();
void createLAMMPSfiles(collagenFibril fib);
#endif
