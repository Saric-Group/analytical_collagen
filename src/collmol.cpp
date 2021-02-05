#include "collmol.hpp"


/* Functions */
void collagenMolecule::readAtoms(std::string &file)
{
  numAtoms = 0;
  std::string line;
  std::ifstream myfile;
  myfile.open(file.c_str(), std::ios::in);
  while (myfile.peek() != EOF) {
    getline(myfile, line);
    numAtoms++;
    charges.push_back(atof(line.c_str()));
    if (charges[numAtoms - 1] != 0) {
      nonZeroChargesIndex.push_back(numAtoms - 1);
    }
  }
  myfile.close();
  length = (numAtoms - 1) * distanceAtoms;
}

void collagenMolecule::readAtoms(std::vector<double> vec)
{
  std::vector<double>().swap(charges);
  std::vector<double>().swap(nonZeroChargesIndex);
  numAtoms = 0;
  for (int i = 0; i < (int) vec.size(); i++) {
    numAtoms++;
    charges.push_back(vec[numAtoms - 1]);
    if (abs(charges[numAtoms - 1]) > 1e-15) {
      nonZeroChargesIndex.push_back(numAtoms - 1);
    }
  }
  length = (numAtoms - 1) * distanceAtoms;
}
