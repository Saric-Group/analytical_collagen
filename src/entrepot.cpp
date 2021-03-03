// // fib.radGap = 34.75;
// // fib.offset = 66.955;
// // cout << "\n\nMinimal energy configuration:";
// // cout << "\nEnergy: " << fib.energy;
// // cout << "\nLateral gap: " << fib.latGap;
// // cout << "\nRadial gap: " << fib.radGap;
// // cout << "\nOffset: " << fib.offset;
// // write3Drescale(fib);
//
// int samples = 500;
// int nativeN = 1054;
// // int nativePos = 86; // * 2 / 4;
// // int nativeNeg = 82 * 2 / 4;
// int nativePos = 22;
// int nativeNeg = 21;
// int pos, neg;
// int error = 0;
// double energy = 0.0;
// double gap = 0.0;
// double offset = 0.0;
// vector<double> dist, dists;
// dists.assign(nativeN, 0.0);
// double targetGap = 34.77;
// double targetP = 66.975;
// double tolerance = 0.1;
// int goodOnes = 0;
// FILE *outf;
// outf = fopen(filePaths.outputpath.c_str(), "a");
// fprintf(outf, "#radGap\toffset\tenergy");
// for (int i = 0; i < samples; i++) {
//   cout << "\nMinimizing " << i + 1 << ". random distribution finished with";
//   dist = createRandomChargeDistribution(nativeN, nativePos, nativeNeg, fib.mol.atomTypes);
//   countCharges(dist, pos, neg);
//   if (pos != nativePos || neg != nativeNeg) {
//     error++;
//   }
//   fib.mol.readAtoms(dist);
//   fib.energy = 1e16;
//   fib.minimizeEnergy();
//   fprintf(outf, "\n%.4f", fib.radGap);
//   fprintf(outf, "\t%.4f", fib.offset);
//   fprintf(outf, "\t%.4f", fib.energy);
//   energy += fib.energy;
//   gap += fib.radGap;
//   offset += fib.offset;
//   cout << " min_energy = " << fib.energy << ".";
//   if (isinf(fib.energy)) {
//     fib.writeXYZ();
//   }
//   if (abs(fib.radGap - targetGap) / targetGap <= tolerance &&
//       abs(fib.offset - targetP) / targetP <= tolerance) {
//     cout << "\t -> within tolerance.";
//     for (int j = 0; j < (int) dist.size(); j++) {
//       dists[j] += dist[j];
//     }
//     goodOnes++;
//   }
// }
//
// for (int j = 0; j < (int) dists.size(); j++) {
//   dists[j] /= goodOnes;
// }
//
// fprintf(outf, "\n\n");
// fprintf(outf, "#Averages:");
// fprintf(outf, "\n%.4f", gap / samples);
// fprintf(outf, "\t%.4f", offset / samples);
// fprintf(outf, "\t%.4f", energy / samples);
// fclose(outf);
// cout << "\n\nErrors: " << error;
//
// file = "./data/random/goodOnes-quarter.dat";
// outf = fopen(file.c_str(), "a");
// fprintf(outf, "#Proportion of distributions with radial gap and");
// fprintf(outf, " offset being within a %.3f tolerance of", tolerance);
// fprintf(outf, " the values for the native distribution:");
// fprintf(outf, " %.3f", 1.0 * goodOnes / samples);
// fprintf(outf, "\n");
// for (int j = 0; j < (int) dists.size(); j++) {
//   fprintf(outf, "%.3f\n", dists[j]);
// }
// fclose(outf);






  // cout << "\nValue: " << filePaths.md_outputpath;
  // cout << "\nValue: " << flags.readCharges;
  // cout << ":: Computation start";

  /************************************************************************/


  // cout << "\n\n:: running computation...";
  /* Only global Emin */
  // fib.mol.readAtoms(filePaths.inputpath);
  // fib.mol.genUniformType();
  // cout << "\nlength=" << fib.mol.length;
  // fib.mol.printAtoms();
  // fib.minimizeEnergy();
  // string file = "./dis/atom_types_1054";
  // fib.mol.readTypes(file);
  // fib.mol.chargesFromTypes();

  // fib.singleEmin();
  // fib.minimizeEnergy();
  // fib.radGap = 34.75;
  // fib.offset = 66.955;
  // cout << "\n\nMinimal energy configuration:";
  // cout << "\nEnergy: " << fib.energy;
  // cout << "\nLateral gap: " << fib.latGap;
  // cout << "\nRadial gap: " << fib.radGap;
  // cout << "\nOffset: " << fib.offset;
  // write3Drescale(fib);
  // genTopologyZero(fib, 6);
  // genInSim(fib);
