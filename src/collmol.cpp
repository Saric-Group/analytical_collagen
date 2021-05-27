#include "collmol.hpp"


/* Functions */
void collagenMolecule::typesFromCharges()
{
  std::vector<int>().swap(atomTypes);
  int type = 0;
  for (int i = 0; i < numAtoms; i++) {
    if (charges[i] < 0) {
      type = maxNumTypes + 1;
    } else if (charges[i] > 0) {
      type = maxNumTypes + 3;
    } else {
      type = maxNumTypes + 2;
    }
    atomTypes.push_back(type);
  }
  numTypes = maxNumTypes + 3;
}
void collagenMolecule::countCharges()
{
  double sum = 0.;
  totalCharge = 0.0;
  numNeg = 0;
  numPos = 0;
  for (int i = 0; i < numAtoms; i++) {
    // std::cout << "\n" << i << "\t" << charges[i];
    sum += charges[i];
    if (charges[i] < 0) numNeg++;
    if (charges[i] > 0) numPos++;
  }
  totalCharge = sum;
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
    if (line.at(0) != '#') {
      numAtoms++;
      charges.push_back(atof(line.c_str()));
      if (charges[numAtoms - 1] != 0) {
        nonZeroChargesIndex.push_back(numAtoms - 1);
      }
    }
  }
  myfile.close();
  length = (numAtoms - 1) * distanceAtoms;
  typesFromCharges();
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
  typesFromCharges();
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
  myfile.close();
  length = (numAtoms - 1) * distanceAtoms;
  chargesFromTypes();
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
  chargesFromTypes();
}

void collagenMolecule::printAtoms()
{
  std::cout << "\n#\n#\n# Atomic information for collagen molecule of length ";
  std::cout << length;
  std::cout << "\n#\tID\tType\tCharge";
  for (int i = 0; i < numAtoms; i++) {
    std::cout << "\n#\t" << i << "\t" << atomTypes[i] << "\t" << charges[i];
  }
}

void collagenMolecule::printMoleculeInfo()
{
  std::cout << "\n#\n#";
  std::cout << "\n# Molecule information:";
  std::cout << "\n#    Number of atoms: " << numAtoms;
  std::cout << "\n#    Number of types: " << numTypes;
  std::cout << "\n#    Molecule length: " << length;
  std::cout << "\n#    Interatomic distance: " << distanceAtoms;
  std::cout << "\n#    Atom diameter: " << diameterAtom;
  std::cout << "\n#    Number of positive charges: " << numPos;
  std::cout << "\n#    Number of negative charges: " << numNeg;
  std::cout << "\n#    Total charge: " << totalCharge;
}

void collagenMolecule::moleculeToFile(std::string &file)
{
  FILE *outf;
  outf = fopen(file.c_str(), "w");
  fprintf(outf, "# charges");
  for (int i = 0; i < numAtoms; i++) {
    fprintf(outf, "\n%.3f", charges[i]);
  }
  fclose(outf);
}

void collagenMolecule::chargesToFile(std::string &file, int k)
{
  FILE *outf;
  outf = fopen(file.c_str(), "a");
  // fprintf(outf, "%i", numAtoms);
  // fprintf(outf, "\nAtoms");
  for (int i = 0; i < numAtoms; i++) {
    fprintf(outf, "\n%.3f", charges[i]);
    fprintf(outf, " %.6f", i * distanceAtoms);
    fprintf(outf, " %.6f", k * diameterAtom);
    fprintf(outf, " 0");
  }
  fclose(outf);
}

void collagenMolecule::addOvitoHeaderToChargeFile(std::string &file)
{
  std::ofstream outputFile("outputFileName");
  std::ifstream inputFile(file);

  std::string line;
  int N = 0;
  while (std::getline(inputFile, line)) {
    N++;
  }
  inputFile.close();

  outputFile << N - 1 << "\n";
  outputFile << "Atoms";

  inputFile.open(file, std::ifstream::in);
  outputFile << inputFile.rdbuf();

  inputFile.close();
  outputFile.close();

  std::remove(file.c_str());
  std::rename("outputFileName", file.c_str());
}

std::vector<double> collagenMolecule::smoothen(int delta)
{
  /* need to redo the boundaries */
  double ma = 0;                  /* moving average */
  std::vector<double> smoothed;        /* smoothed moving average charges */
  double div = 2.0 * delta + 1;

  // std::cout << "######## Start ########" << std::endl;
  for (int i = -delta; i <= delta; i++) {
    if (i < 0) {
      ma += charges[abs(i) - 1] / div;
      // std::cout << "+" << abs(i) - 1;
      // std::cout << " +" << charges[abs(i) - 1] / div << std::endl;
    } else {
      ma += charges[i] / div;
      // std::cout << "+" << i;
      // std::cout << " +" << charges[i] / div << std::endl;
    }
  }
  // std::cout << "moving average for site 0 " << ma << std::endl;
  smoothed.push_back(ma);
  // std::cout << "######## Body ########" << std::endl;
  for (int i = 1; i < numAtoms; i++) {
    // std::cout << "i = " << i << std::endl;
    if (i - 1 - delta < 0) {
      ma += charges[i + delta] / div;
      // std::cout << "+" << i + delta;
      // std::cout << " +" << charges[i + delta] / div;
      ma -= charges[abs(i - delta)] / div;
      // std::cout << " -" << abs(i - delta);
      // std::cout << " -" << charges[abs(i - delta)] / div << std::endl;
    } else if (i + delta > numAtoms - 1) {
      ma -= charges[i - 1 - delta] / div;
      // std::cout << "-" << i - 1 - delta;
      // std::cout << " -" << charges[i - 1 - delta] / div;
      ma += charges[numAtoms - 1 - (i - numAtoms + delta)] / div;
      // std::cout << " +" << numAtoms - 1 - (i - numAtoms + delta);
      // std::cout << " +" << charges[numAtoms - 1 - (i - numAtoms + delta)] / div << std::endl;
    } else {
      ma += (charges[i + delta] - charges[i - 1 - delta]) / div;
      // std::cout << "+" << i + delta;
      // std::cout << " +" << charges[i + delta] / div;
      // std::cout << " -" << i - 1 - delta;
      // std::cout << " -" << charges[i - 1 - delta] / div << std::endl << std::endl;
    }
    // std::cout << "moving average for site " << i << " " << ma << std::endl;
    smoothed.push_back(ma);
    // if (i >= 12050) {
    //   std::cin >> i;
    // }
  }

  // printVec(smoothed);
  // double sum = 0.0;
  // for (int i = 0; i < (int) smoothed.size(); i++) {
  //   sum += smoothed[i];
  // }
  // std::cout << "\nSize: " << smoothed.size();
  // std::cout << "\nSum after smoothing: " << sum;

  return smoothed;
}

std::vector<double> collagenMolecule::binNormalize(std::vector<double> vec, int n)
{
  std::vector<double> normalized;
  double ma;
  int binSize = floor(numAtoms / n);
  int split = numAtoms % n;
  for (int i = 0; i < split; i++) {
    // cout << "\nBin " << i + 1 << " adding:\n";
    ma = 0;
    for (int j = 0; j <= binSize; j++) {
      // cout << i * (binSize + 1) + j << " ";
      ma += vec[i * (binSize + 1) + j];
    }
    normalized.push_back(ma);
  }
  for (int i = split; i < n - 1; i++) {
    // cout << "\nBin " << i + 1 << " adding:\n";
    ma = 0;
    for (int j = 0; j < binSize; j++) {
      // cout << i * binSize + split + j << " ";
      ma += vec[i * binSize + split + j];
    }
    normalized.push_back(ma);
  }
  ma = 0;
  // cout << "\nBin " << n << " adding:\n";
  for (int i = (n - 1) * binSize + split; i < (int) vec.size(); i++) {
    // cout << i << " ";
    ma += vec[i];
  }
  normalized.push_back(ma);

  // cout << "\n\ntest: " << vec.size() % n;

  // double sum = 0.0;
  // for (int i = 0; i < (int) normalized.size(); i++) {
  //   sum += normalized[i];
  // }
  // cout << "\nSize: " << normalized.size();
  // std::cout << "\nIntegral: " << sum;

  return normalized;
}
