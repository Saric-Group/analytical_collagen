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
    double x, y, z;
    position pos;
    for (int i = 0; i < numPoints; i++) {
      x = i * length / numPoints - length / 2.0;
      for (int j = 0; j < numPoints; j++) {
        y = j * length / numPoints - length / 2.0;
        for (int k = 0; k < numPoints; k++) {
          z = k * length / numPoints - length / 2.0;
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
#endif
