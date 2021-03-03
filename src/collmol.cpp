#include "collmol.hpp"


/* Functions */
void collagenMolecule::countCharges()
{
  double sum = 0.;
  for (int i = 0; i < numAtoms; i++) {
    // std::cout << "\n" << i << "\t" << charges[i];
    sum += charges[i];
    if (charges[i] < 0) numNeg++;
    if (charges[i] > 0) numPos++;
  }
  totalCharge = sum;
}
void collagenMolecule::readCharges(std::string &file)
{
  std::vector<double>().swap(charges);
  std::vector<double>().swap(nonZeroChargesIndex);
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
  countCharges();
}

void collagenMolecule::readCharges(std::vector<double> vec)
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
  countCharges();
}

void collagenMolecule::readTypes(std::string &file)
{
  std::vector<int>().swap(atomTypes);
  numAtoms = 0;
  std::string line;
  std::ifstream myfile;
  myfile.open(file.c_str(), std::ios::in);
  while (myfile.peek() != EOF) {
    getline(myfile, line);
    if (line.at(0) != '#') {
      numAtoms++;
      atomTypes.push_back(atoi(line.c_str()));
    }
  }
  length = (numAtoms - 1) * distanceAtoms;
}

void collagenMolecule::readTypes(std::vector<int> vec)
{
  std::vector<int>().swap(atomTypes);
  numAtoms = 0;
  for (int i = 0; i < (int) vec.size(); i++) {
    atomTypes.push_back(vec[i]);
    numAtoms++;
  }
  length = (numAtoms - 1) * distanceAtoms;
}

void collagenMolecule::chargesFromTypes()
{
  std::vector<double>().swap(charges);
  std::vector<double>().swap(nonZeroChargesIndex);
  // numAtoms = atomTypes.size();
  for (int i = 0; i < (int) atomTypes.size(); i++) {
    // std::cout << "\n" << i << " " << typeCharge[atomTypes[i] - 1];
    charges.push_back(chargeAtom * typeCharge[atomTypes[i] - 1]);
    if(abs(charges[i]) > 1e-15) {
      nonZeroChargesIndex.push_back(i);
    }
  }
  countCharges();
}

void collagenMolecule::genUniformType()
{
  std::vector<int>().swap(atomTypes);
  for (int i = 0; i < numAtoms; i++) {
    atomTypes.push_back(1);
  }
  numTypes = 1;
}

void collagenMolecule::printAtoms()
{
  std::cout << "\n#\n# Atomic information for collagen molecule of length ";
  std::cout << length;
  std::cout << "\n#\tID\tType\tCharge";
  for (int i = 0; i < numAtoms; i++) {
    std::cout << "\n#\t" << i << "\t" << atomTypes[i] << "\t" << charges[i];
  }
}
