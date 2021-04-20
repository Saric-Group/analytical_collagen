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

void runRandomAnalysis(int samples, layerModel lm)
{
  // double targetGap = 34.77;
  double targetDper = 66.975;
  double tolerance = 0.15;

  double binSize = 2.5;
  int nBins = ceil(lm.mol.length / binSize);
  std::vector<int> binsGap, binsGoodGap, binsDper;
  std::vector<std::vector<int>> bins2D;
  binsGap.assign(nBins, 0);
  binsGoodGap.assign(nBins, 0);
  binsDper.assign(nBins, 0);
  bins2D.assign(nBins, std::vector<int> (nBins, 0));

  FILE *outf;
  std::string filepath = filePaths.outputpath;
  filepath += "_N=" + std::to_string(lm.mol.numAtoms);
  filepath += "_nPos=" + std::to_string(lm.mol.numPos);
  filepath += "_nNeg=" + std::to_string(lm.mol.numNeg);
  filepath += "_layers=" + std::to_string(lm.layers);
  outf = fopen((filepath + ".dat").c_str(), "a");
  fprintf(outf, "#radGap\toffset\tenergy");

  std::cout << "\n#";
  for (int i = 0; i < samples; i++) {
    std::cout << "\n# Sample " << i + 1 << " / " << samples;
    std::cout.flush();
    lm.mol.readCharges(createRandomChargeDistribution(lm.mol));
    lm.energy = 1e16;
    lm.minimizeEnergy();
    fprintf(outf, "\n%.4f", lm.radGap);
    fprintf(outf, "\t%.4f", lm.offset);
    fprintf(outf, "\t%.4f", lm.energy);
    binsGap[(int) lm.radGap / binSize]++;
    binsDper[(int) lm.offset / binSize]++;
    bins2D[(int) lm.radGap / binSize][(int) lm.offset / binSize]++;
    if (abs(lm.offset - targetDper) / targetDper <= tolerance) {
      binsGoodGap[(int) lm.radGap / binSize]++;
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

  outf = fopen((filepath + ".2Dbins").c_str(), "w");
  fprintf(outf, "#gap\tdper\tcounts");
  for (int i = 0; i < nBins; i++) {
    for (int j = 0; j < nBins; j++) {
      fprintf(outf, "\n%.3f", binSize * i);
      fprintf(outf, "\t%.3f", binSize * j);
      fprintf(outf, "\t%i", bins2D[i][j]);
    }
    fprintf(outf, "\n");
  }
  fclose(outf);

  outf = fopen((filePaths.outputpath + ".peaks").c_str(), "a");
  fprintf(outf, "\n");
  fprintf(outf, "%i", lm.mol.numPos);
  fprintf(outf, "\t%i", lm.mol.numNeg);
  fprintf(outf, "\t%i", lm.mol.numAtoms);
  fprintf(outf, "\t%.1f", binSize * index_of_max_dper);
  fprintf(outf, "\t%i", binsDper[index_of_max_dper]);
  fprintf(outf, "\t%.1f", binSize * index_of_max_gap);
  fprintf(outf, "\t%i", binsGap[index_of_max_gap]);
  fprintf(outf, "\t%i", samples);
  fclose(outf);
}
