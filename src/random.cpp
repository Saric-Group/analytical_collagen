#include "random.hpp"


/* Variables */
// int seed = 1337;
std::random_device device;
std::mt19937 generator(device());
// std::mt19937 generator(seed);
extern filePaths_ filePaths;


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

std::vector<double> createRandomChargeDistribution(collagenMolecule mol)
{
  std::uniform_int_distribution<int> dis_int(0, mol.numAtoms - 1);
  std::vector<double> vec;
  vec.assign(mol.numAtoms, 0);
  int pos;
  int counter = 0;
  bool exit;

  while (counter < std::max(mol.numPos, mol.numNeg)) {
    if (counter < mol.numPos) {
      exit = false;
      while (!exit) {
        pos = dis_int(generator);
        if (vec[pos] == 0 && mol.atomTypes[pos] != 5) {
          vec[pos] += 22.4;
          exit = true;
        }
      }
    }
    if (counter < mol.numNeg) {
      exit = false;
      while (!exit) {
        pos = dis_int(generator);
        if (vec[pos] == 0 && mol.atomTypes[pos] != 5) {
          vec[pos] -= 22.4;
          exit = true;
        }
      }
    }
    counter++;
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

void runRandomAnalysis(int samples, collagenFibril fib)
{
  // double targetGap = 34.77;
  double targetDper = 66.975;
  double tolerance = 0.15;

  double binSize = 2.5;
  int nBins = ceil(fib.mol.length / binSize);
  std::vector<int> binsGap, binsGoodGap, binsDper;
  binsGap.assign(nBins, 0);
  binsGoodGap.assign(nBins, 0);
  binsDper.assign(nBins, 0);

  FILE *outf;
  std::string filepath = filePaths.outputpath;
  filepath += "_nPos=" + std::to_string(fib.mol.numPos);
  filepath += "_nNeg=" + std::to_string(fib.mol.numNeg);
  outf = fopen((filepath + ".dat").c_str(), "w");
  fprintf(outf, "#radGap\toffset\tenergy");

  for (int i = 0; i < samples; i++) {
    std::cout << "\n# Run " << i + 1;
    std::cout.flush();
    fib.mol.readCharges(createRandomChargeDistribution(fib.mol));
    fib.energy = 1e16;
    fib.minimizeEnergy();
    fprintf(outf, "\n%.4f", fib.radGap);
    fprintf(outf, "\t%.4f", fib.offset);
    fprintf(outf, "\t%.4f", fib.energy);
    binsGap[(int) fib.radGap / binSize]++;
    binsDper[(int) fib.offset / binSize]++;
    if (abs(fib.offset - targetDper) / targetDper <= tolerance) {
      binsGoodGap[(int) fib.radGap / binSize]++;
    }
  }
  fclose(outf);

  int index_of_max_dper = 0, index_of_max_gap = 0;
  outf = fopen((filepath + ".bins").c_str(), "w");
  fprintf(outf, "#bin\tdper\tgap\tgoodGap");
  for (int i = 0; i < nBins; i++) {
    fprintf(outf, "\n%.3f", binSize * i);
    fprintf(outf, "\t%i", binsDper[i]);
    fprintf(outf, "\t%i", binsGap[i]);
    fprintf(outf, "\t%i", binsGoodGap[i]);
    if (binsDper[i] > binsDper[index_of_max_dper]) {
      index_of_max_dper = i;
    }
    if (binsGap[i] > binsGap[index_of_max_gap]) {
      index_of_max_gap = i;
    }
  }
  fclose(outf);

  outf = fopen((filePaths.outputpath + ".peaks").c_str(), "a");
  fprintf(outf, "\n");
  fprintf(outf, "%i", fib.mol.numPos);
  fprintf(outf, "\t%i", fib.mol.numNeg);
  fprintf(outf, "\t%i", fib.mol.numAtoms);
  fprintf(outf, "\t%.1f", binSize * index_of_max_dper);
  fprintf(outf, "\t%i", binsDper[index_of_max_dper]);
  fprintf(outf, "\t%.1f", binSize * index_of_max_gap);
  fprintf(outf, "\t%i", binsGap[index_of_max_gap]);
  fprintf(outf, "\t%i", samples);
  fclose(outf);
}
