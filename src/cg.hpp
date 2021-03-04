#ifndef CG
#define CG

#include "main.hpp"


/* Functions */
std::vector<double> smoothenSMCA(int delta,
                                 std::vector<double> vec,
                                 bool avg = false);
std::vector<double> binNormalize(int n, std::vector<double> vec);

#endif
