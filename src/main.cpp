#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"
#include "random.hpp"

using namespace std;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;


/* Settings */
filePaths_ filePaths;       /* IO File paths */
flags_ flags;               /* Various flags, see main.hpp */

int main(int argc, char const *argv[])
{

  collagenFibril fib;

  if (parse_all_args(argc, argv, fib) != 1) {
    return 0;
  }

  time_point start = chrono::high_resolution_clock::now();
  cout << ":: Computation start";

  /************************************************************************/

  cout << "\n\n:: running computation...";
  /* Only global Emin */
  fib.mol.readAtoms(filePaths.inputpath);
  fib.singleEmin();

  /************************************************************************/

  time_point end = chrono::high_resolution_clock::now();
  chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
  cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";

  /************************************************************************/

  return 0;
}


/* Functions */
// void readTypes(string &file)
// {
//   int counter = 0;
//   string line;
//   ifstream myfile;
//   myfile.open(file.c_str(), ios::in);
//   while (myfile.peek() != EOF) {
//     getline(myfile, line);
//     if (line.at(0) != '#') {
//       counter++;
//       type.push_back(atoi(line.c_str()));
//       if (type[counter - 1] == 2 ||
//           type[counter - 1] == 10 ||
//           type[counter - 1] == 13 ||
//           type[counter - 1] == 15 ||
//           type[counter - 1] == 20) {
//           type_ind.push_back(counter - 1);
//       }
//     }
//   }
//   cout << "\n -> " << type_ind.size() << " hydrophobic atoms found.";
// }

// double newLJ_per_mol(double pos, double dx, double ref, double lat_gap)
// {
//   double d, sum;
//   int left, right, atomtype;
//   sum = 0;
//
//   left = ceil(((pos - dx) - ref) / distance_atoms);
//   right = floor(((pos + dx) - ref) / distance_atoms);
//
//   if (left < N && right >= 0) {
//       left = max(left, 0);
//       right = min(right, N - 1);
//       for (int i = 0; i <= right - left; i++) {
//           atomtype = type[left + i];
//           if (atomtype == 2 ||
//               atomtype == 10 ||
//               atomtype == 13 ||
//               atomtype == 15 ||
//               atomtype == 20) {
//               d = distance(pos, ref, left + i, lat_gap);
//               // cout << "\nhit with " << factorCD(q1, q2, d) / 150;
//               // cout << " to change sum = " << sum / 150;
//               sum += factorLJ(d);
//               // cout << " to sum = " << sum / 150;
//           }
//       }
//   }
//   return sum;
// }

// double newLJ_layer(double pos, double dx, int layer, double offset,
//                 double lat_gap)
// {
//   double sum = 0;
//   sum += newLJ_per_mol(pos, dx, layer * offset - box, lat_gap);
//   sum += newLJ_per_mol(pos, dx, layer * offset, lat_gap);
//   sum += newLJ_per_mol(pos, dx, layer * offset + box, lat_gap);
//
//   return sum;
// }

// double newLJ(double lat_gap, double rad_gap, double offset, int layers)
// {
//   int atomtype, left, right;
//   double pos, sum, d;
//   double dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
//   sum = 0;
//   box = L + rad_gap;
//   for (int atom = 0; atom < (int) type_ind.size(); atom++) {
//     pos = type_ind[atom] * distance_atoms;
//
//     /* Atom atom feels from molecules from layer i: */
//     for (int layer = 1; layer < layers; layer++) {
//       sum += newLJ_layer(pos, dx, layer, offset, lat_gap);
//     }
//
//     /* Atom atom feels from left molecule: */
//     left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
//     right = floor(((pos + lj_cutoff) + box) / distance_atoms);
//     if (left < N && right >= 0) {
//         left = max(left, 0);
//         right = min(right, N - 1);
//         for (int i = 0; i <= right - left; i++) {
//             atomtype = type[left + i];
//             if (atomtype == 2 ||
//                 atomtype == 10 ||
//                 atomtype == 13 ||
//                 atomtype == 15 ||
//                 atomtype == 20) {
//                 d = abs(pos + box - (left + i) * distance_atoms);
//                 sum += factorLJ(d);
//             }
//         }
//     }
//
//     /* Atom atom feels from right molecule: */
//     left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
//     right = floor(((pos + cd_cutoff) - box) / distance_atoms);
//     if (left < N && right >= 0) {
//         left = max(left, 0);
//         right = min(right, N - 1);
//         for (int i = 0; i <= right - left; i++) {
//             atomtype = type[left + i];
//             if (atomtype == 2 ||
//                 atomtype == 10 ||
//                 atomtype == 13 ||
//                 atomtype == 15 ||
//                 atomtype == 20) {
//                 d = abs(pos - box - (left + i) * distance_atoms);
//                 sum += factorLJ(d);
//             }
//         }
//     }
//   }
//   return sum;
// }


/* this function is not maintained anymore */
// void multipleEmin(string &file, int layers, int number)
// {
//   double lat_gap = diameter_atom;
//   double rad_gap, offset;
//
//   int rad_max = (int) ceil(L / (layers - 1.0) / distance_atoms);
//   int off_max;
//   double lj, lje;
//   double cd, cde;
//   vector<vector<vector<double> > > latmin, radmin, offmin;
//   vector<vector<vector<double> > > eCDmin, eLJmin, eTotmin;
//
//   /* Write header for outputfile */
//   string tmps = file;
//   for (int i = 1; i <= number; i++) {
//     file = tmps + "Emin_" + to_string(i) + ".dat";
//     header(file);
//   }
//
//   /* Assign needed data vectors */
//   // latmin.assign(10, vector<vector<double> > (50, vector<double> (20, 0)));
//   radmin.assign(number, vector<vector<double> > (lj_steps, vector<double> (cd_steps, 0)));
//   offmin.assign(number, vector<vector<double> > (lj_steps, vector<double> (cd_steps, 0)));
//   eLJmin.assign(number, vector<vector<double> > (lj_steps, vector<double> (cd_steps, 1e6)));
//   eCDmin.assign(number, vector<vector<double> > (lj_steps, vector<double> (cd_steps, 1e6)));
//   eTotmin.assign(number, vector<vector<double> > (lj_steps,
//                                                   vector<double> (cd_steps, 1e6)));
//
//   for (int rad = 0; rad <= rad_max; rad++) {
//     cout << "\n --> " << (100. * rad) / (1. * rad_max);
//     cout << "\% done ...";
//     rad_gap = rad * distance_atoms;
//     double tmp = (L - (layers - 1.0) * rad_gap) / layers / distance_atoms;
//     off_max = (int) ceil(tmp);
//     for (int off = 0; off <= off_max; off++) {
//       offset = rad_gap + off * distance_atoms;
//       /* Energy factors calculation for given parameters above */
//       lj = compLJ(lat_gap, rad_gap, offset, layers);
//       cd = compCD(lat_gap, rad_gap, offset, layers);
//       /* Both CD and LJ potential */
//       for (int i = 0; i < lj_steps; i++) {
//           lje = (lj_min + (i) * lj_stepsize) * lj / 5;
//           for (int j = 0; j < cd_steps; j++) {
//               cde = cd / (cd_min + (j) * cd_stepsize / 2);
//
//               for (int k = 0; k < number; k++) {
//                 if (lje + cde < eTotmin[k][i][j]) {
//                     eLJmin[k][i][j] = lje;
//                     eCDmin[k][i][j] = cde;
//                     eTotmin[k][i][j] = lje + cde;
//                     radmin[k][i][j] = rad_gap;
//                     offmin[k][i][j] = offset;
//                     break;
//                 }
//               }
//           }
//       }
//     }
//   }
//
//   FILE *outf;
//   for (int l = 1; l <= number; l++) {
//     int k = l - 1;
//     file = tmps + "Emin_" + to_string(l) + ".dat";
//     outf = fopen(file.c_str(), "a");
//     /* Both LJ and CD potential */
//     for (int i = 0; i < lj_steps; i++) {
//         for (int j = 0; j < cd_steps; j++) {
//             fprintf(outf, "\n");
//             fprintf(outf, "%.3f", lj_min + (i) * lj_stepsize / 5.);
//             fprintf(outf, "\t%.3f", cd_min + (j) * cd_stepsize / 2.);
//             fprintf(outf, "\t%.3f", lat_gap);
//             fprintf(outf, "\t%.3f", radmin[k][i][j]);
//             fprintf(outf, "\t%.3f", offmin[k][i][j]);
//             fprintf(outf, "\t%.3f", L + radmin[k][i][j] - offmin[k][i][j]);
//             fprintf(outf, "\t%.3f", eCDmin[k][i][j]);
//             fprintf(outf, "\t%.3f", eLJmin[k][i][j]);
//             fprintf(outf, "\t%.3f", eTotmin[k][i][j]);
//         }
//         fprintf(outf, "\n");
//     }
//   }
//   fclose(outf);
// }

// double energy(double lat_gap, double rad_gap, double offset, int layer1,
//               int layer2, double lj, double cd, bool lj_on)
// {
//   if (layer1 >= layer2) {
//     cout << "error";
//     return 0.0;
//   }
//
//   double pos, sum, q1;
//   double dx = sqrt(cd_cutoff * cd_cutoff - lat_gap * lat_gap);
//   sum = 0.0;
//   box = L + rad_gap;
//   /* cycle through atoms of layer1 */
//   for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
//     q1 = charge[charge_ind[atom]];
//     pos = charge_ind[atom] * distance_atoms;
//     /* focus atom1 feels from layer2 */
//     sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset - box, lat_gap);
//     sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset, lat_gap);
//     sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset + box, lat_gap);
//   }
//   sum *= cd;
//
//   /* turn on/off for hydrophobic interactions */
//   if (lj_on) {
//     dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
//     for (int atom = 1; atom <= N; atom++) {
//       pos = (atom - 1) * distance_atoms;
//       sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset - box, lat_gap);
//       sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset, lat_gap);
//       sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset + box, lat_gap);
//     }
//   }
//   return sum;
// }
//
// double intra_layer_energy(double lat_gap, double rad_gap, double offset,
//                           double lj, double cd, bool lj_on)
// {
//   double pos, sum, q1, q2, d;
//   int left, right;
//   sum = 0;
//   box = L + rad_gap;
//   for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
//     q1 = charge[charge_ind[atom]];
//     pos = charge_ind[atom] * distance_atoms;
//
//     /* Atom atom feels from left molecule: */
//     left = ceil(((pos - cd_cutoff) + box) / distance_atoms);
//     right = floor(((pos + cd_cutoff) + box) / distance_atoms);
//     if (left < N && right >= 0) {
//         left = max(left, 0);
//         right = min(right, N - 1);
//         for (int i = 0; i <= right - left; i++) {
//             q2 = charge[left + i];
//             if (q2 != 0) {
//                 d = abs(pos + box - (left + i) * distance_atoms);
//                 sum += factorCD(q1, q2, d);
//             }
//         }
//     }
//
//     /* Atom atom feels from right molecule: */
//     left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
//     right = floor(((pos + cd_cutoff) - box) / distance_atoms);
//     if (left < N && right >= 0) {
//         left = max(left, 0);
//         right = min(right, N - 1);
//         for (int i = 0; i <= right - left; i++) {
//             q2 = charge[left + i];
//             if (q2 != 0) {
//                 d = abs(pos - box - (left + i) * distance_atoms);
//                 sum += factorCD(q1, q2, d);
//             }
//         }
//     }
//   }
//   sum *= cd;
//
//   if (lj_on) {
//     for (int atom = 1; atom <= N; atom++) {
//       /* Atom atom feels from left molecule: */
//       left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
//       right = floor(((pos + lj_cutoff) + box) / distance_atoms);
//       if (left < N && right >= 0) {
//           left = max(left, 0);
//           right = min(right, N - 1);
//           for (int i = 0; i <= right - left; i++) {
//               double d = abs(pos + box - (left + i) * distance_atoms);
//               sum += lj * factorLJ(d);
//           }
//       }
//
//       /* Atom atom feels from right molecule: */
//       left = ceil(((pos - lj_cutoff) - box) / distance_atoms);
//       right = floor(((pos + lj_cutoff) - box) / distance_atoms);
//       if (left < N && right >= 0) {
//           left = max(left, 0);
//           right = min(right, N - 1);
//           for (int i = 0; i <= right - left; i++) {
//               double d = abs(pos - box - (left + i) * distance_atoms);
//               sum += lj * factorLJ(d);
//           }
//       }
//     }
//   }
//
//   return sum;
// }
