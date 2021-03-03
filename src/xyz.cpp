#include "xyz.hpp"


/* Variables */
extern filePaths_ filePaths;

/* Functions */
void write3Drescale(collagenFibril fib)
{
  double res = 10.0;
  double box = fib.radGap + fib.mol.length;

  FILE *outf;
  std::string file = filePaths.outputpath;
  if (file.find(filePaths.file_extension) != std::string::npos)
  {
      file = file.substr(0, file.find(filePaths.file_extension));
  }
  file += ".xyz";
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "%i", 1 * 3 * fib.mol.numAtoms);
  fprintf(outf, "\n");

  /* focus */
  for (int i = 0; i < fib.mol.numAtoms; i++) {
    fprintf(outf, "\n0");
    fprintf(outf, "\t%.8f", i * fib.mol.distanceAtoms / res);
    fprintf(outf, "\t%.8f", 0.);
    fprintf(outf, "\t%.8f", 0.);
  }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  // }
  // //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  // }
  //
  //
  //
  //
  for (int i = 0; i < fib.mol.numAtoms; i++) {
    fprintf(outf, "\n5");
    fprintf(outf, "\t%.8f", (box + i * fib.mol.distanceAtoms) / res);
    fprintf(outf, "\t%.8f", 0.);
    fprintf(outf, "\t%.8f", 0.);
  }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  // }
  // //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  // }
  //
  //
  //
  //
  //
  for (int i = 0; i < fib.mol.numAtoms; i++) {
    fprintf(outf, "\n5");
    fprintf(outf, "\t%.8f", (-box + i * fib.mol.distanceAtoms) / res);
    fprintf(outf, "\t%.8f", 0.);
    fprintf(outf, "\t%.8f", 0.);
  }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (-box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n2");
  //   fprintf(outf, "\t%.8f", (-box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  //   fprintf(outf, "\t%.8f", 0.);
  // }
  //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (-box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", fib.mol.diameterAtom);
  // }
  // //
  // for (int i = 0; i < fib.mol.numAtoms; i++) {
  //   fprintf(outf, "\n4");
  //   fprintf(outf, "\t%.8f", (-box + 0 * fib.offset + i * fib.mol.distanceAtoms) / res);
  //   fprintf(outf, "\t%.8f", 0.);
  //   fprintf(outf, "\t%.8f", -fib.mol.diameterAtom);
  // }

  fclose(outf);
}
