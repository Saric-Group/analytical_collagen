#ifndef RANDOM
#define RANDOM

#include "main.hpp"
#include "layermodel.hpp"

#include <random>


/* Functions */
std::vector<double> createRandomChargeDistribution(int n, int nPos, int nNeg, std::vector<int> types);
std::vector<double> createRandomChargeDistribution(collagenMolecule mol);
void countCharges(std::vector<double> vec, int &pos, int&neg);
void printChargeDistribution(std::vector<double> vec);

void runRandomAnalysis(int samples, layerModel lm);


#endif
