#include "random.hpp"


/* Variables */
// int seed = 1337;
std::random_device device;
std::mt19937 generator(device());
// std::mt19937 generator(seed);


/* Functions */
std::vector<double> createRandomChargeDistribution(int n, int nPos, int nNeg, std::vector<int> types)
{
  std::uniform_int_distribution<int> dis_int(0, n - 1);
  std::vector<double> vec;
  vec.assign(n, 0);
  int pos;
  int counter = 0;

  while (counter < nPos) {
    pos = dis_int(generator);
    if (vec[pos] == 0 && types[pos] != 5) {
      vec[pos] += 22.4;
      counter++;
    }
  }
  counter = 0;
  while (counter < nNeg) {
    pos = dis_int(generator);
    if (vec[pos] == 0 && types[pos] != 5) {
      vec[pos] -= 22.4;
      counter++;
    }
  }

  return vec;
}

void countCharges(std::vector<double> vec, int &pos, int&neg)
{
  pos = 0;
  neg = 0;
  for (int i = 0; i < (int) vec.size(); i++) {
    if (vec[i] > 0) pos++;
    if (vec[i] < 0) neg++;
  }
}

void printChargeDistribution(std::vector<double> vec)
{
  for (int i = 0; i < (int) vec.size(); i++) {
    std::cout << vec[i] << " ";
  }
}
