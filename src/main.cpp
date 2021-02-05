#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"

using namespace std;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;

/* Settings */
parameters_ parameters;   /* Calculation parameters */
filePaths_ filePaths; /* File settings */
flags_ flags;             /* Flags */



// int produce_xyz(int numatoms,int rows,double distatom, double latgap, double radgap, double offset ,int id[],int type[],double charges[], double xpos[], double ypos[], double zpos[]){
//   int natoms = layers*rows*numatoms;
//
//   //cout << natoms << endl;
//   //cout << "Atoms. Timestep: 0" << endl;
//
//   map<int,int> typemap;
//   int typecount = 0;
//   for(int i=0; i<numatoms; i++){
//     if(not typemap.count(charge[i])){
//       typemap[charge[i]]=typecount;
//       typecount += 1;
//     }
//   }
//
//   for(int l=0; l<layers; l++){
//     for(int r=0; r<rows; r++){
//       for(int a=0; a<numatoms; a++){
//         int i = l*rows*numatoms + r*numatoms + a;
//         double drows = (double) rows;
//         double dr =(double) r;
//         double da =(double) a;
//         double dl =(double) l;
//         id[i] = i;
//         charges[i] = charge[a];
//         type[i] = typemap[charge[a]];
//         xpos[i] = (((double)numatoms)*dr*distatom) + (dr*radgap) + dl*offset + da*distatom;
//         ypos[i] = -dl*latgap;
//         zpos[i] = 0.0;
//         //cout << id[i] << "  " << xpos[i] << "  " << ypos[i] << "  " << zpos[i] << "  " << type[i] << "  " << charges[i] << endl;
//       }
//     }
//   }
//
//   return 1;
//
// }
//
// int output_xyz(string filename, int rows, double latgap, double radgap, double offset){
//
//   int natoms = layers*rows*N;
//   int * id = new int[natoms];
//   int * type = new int[natoms];
//   double * charges = new double[natoms];
//   double * xpos = new double[natoms];
//   double * ypos = new double[natoms];
//   double * zpos = new double[natoms];
//
//   int err = produce_xyz(N,rows,distance_atoms,latgap,radgap,offset,id,type,charges,xpos,ypos,zpos);
//
//   if(err != 1){
//     return 0;
//   }
//
//   FILE *xyzf;
//   xyzf = fopen(filename.c_str(), "w");
//
//   fprintf(xyzf, "%i\n", natoms);
//   fprintf(xyzf, "Atoms. Timestep: 0\n");
//   for(int i=0; i<natoms;i++){
//     string xyzrow =
//     to_string(id[i])
//     + "  "
//     + clean_to_string(xpos[i])
//     + "  "
//     + clean_to_string(ypos[i])
//     + "  "
//     + clean_to_string(zpos[i])
//     + "  "
//     + to_string(type[i])
//     + "  "
//     + clean_to_string(charges[i])
//     + "\n";
//     fprintf(xyzf, "%s",xyzrow.c_str());
//   }
//   fclose (xyzf);
//
//   return err;
//
// }



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
  cout << "\n layers: " << fib.layers;
  cout << "\n distanceAtoms: " << fib.mol.distanceAtoms;
  cout << "\n atoms: " << fib.mol.numAtoms << "\n";
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

void header(string &file)
{
  FILE *outf;
  outf = fopen(file.c_str(), "w");
  fprintf(outf, "#atoms per molecule N = %i", N);
  fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
  fprintf(outf, "\n#molecule_length = %.3f", L);
  fprintf(outf, "\n#layers = %i", layers);
  fprintf(outf, "\n\n#LJ_epsilon");
  fprintf(outf, "\tCD_epsilon");
  fprintf(outf, "\tlateral_gap_min");
  fprintf(outf, "\tradial_gap_min");
  fprintf(outf, "\toffset_min");
  fprintf(outf, "\tD-periodicity_min");
  fprintf(outf, "\tE_CD_min");
  fprintf(outf, "\tE_LJ_min");
  fprintf(outf, "\tE_total_min");
  fclose(outf);
}

// void singleEmin(string &file, int layers, collagenFibril fib)
// {
//   double lat_gap = diameter_atom;
//   double rad_gap, offset;
//
//   // int rad_max = (int) ceil(L / (layers - 1.0) / distance_atoms);
//   int rad_max = (int) ceil((L + 2. * max_cutoff) / 0.1);
//   int off_max;
//
//   double lj, lje;
//   double cd, cde;
//   vector<vector<double> > latmin, radmin, offmin, eCDmin, eLJmin, eTotmin;
//
//
//
//   /* Assign needed data vectors */
//   latmin.assign(lj_steps, vector<double> (cd_steps, 0));
//   radmin.assign(lj_steps, vector<double> (cd_steps, 0));
//   offmin.assign(lj_steps, vector<double> (cd_steps, 0));
//   eLJmin.assign(lj_steps, vector<double> (cd_steps, 1e6));
//   eCDmin.assign(lj_steps, vector<double> (cd_steps, 1e6));
//   eTotmin.assign(lj_steps, vector<double> (cd_steps, 1e6));
//
//   int lat_max = (int) ceil(max_cutoff / (0.1 * diameter_atom));
//   for (int lat = 0; lat <= lat_max; lat++) {
//     lat_gap = lat * 0.1 * diameter_atom;
//   for (int rad = 0; rad <= rad_max; rad++) {
//
//     rad_gap = rad * distance_atoms;
//     double tmp = (L - (layers - 1.0) * rad_gap) / layers / distance_atoms;
//     off_max = (int) ceil(tmp);
//
//     for (int off = 0; off <= off_max; off++) {
//       offset = rad_gap + off * distance_atoms;
//       /* Energy factors calculation for given parameters above */
//       lj = compLJ(lat_gap, rad_gap, offset, layers);
//       cd = compCD(lat_gap, rad_gap, offset, layers);
//       /* Both CD and LJ potential */
//       for (int i = 0; i < lj_steps; i++) {
//           lje = (lj_min + (i) * lj_stepsize) * lj;
//           for (int j = 0; j < cd_steps; j++) {
//               cde = cd / (cd_min + (j) * cd_stepsize);
//               if (lje + cde < eTotmin[i][j]) {
//                   eLJmin[i][j] = lje;
//                   eCDmin[i][j] = cde;
//                   eTotmin[i][j] = lje + cde;
//                   latmin[i][j] = lat_gap;
//                   radmin[i][j] = rad_gap;
//                   offmin[i][j] = offset;
//               }
//           }
//       }
//     }
//   }
//   }
//
//
//   if(charge_hashed_outputs)
//   {
//     hash_output(fib);
//   }
//
//   if (file.find(file_extension) != string::npos)
//   {
//       file = file.substr(0, file.find(file_extension));
//   }
//   string xyzbasefile = "";
//
//   if(xyz_outputs){
//     xyzbasefile = file;
//   }
//
//   file+=".dat";
//
//   /* Write header for outputfile */
//   header(file);
//
//
//   FILE *outf;
//   outf = fopen(file.c_str(), "a");
//   /* Both LJ and CD potential */
//   for (int i = 0; i < lj_steps; i++) {
//       for (int j = 0; j < cd_steps; j++) {
//           fprintf(outf, "\n");
//           fprintf(outf, "%.3f", lj_min + (i) * lj_stepsize);            //LJ_epsilon
//           fprintf(outf, "\t%.3f", cd_min + (j) * cd_stepsize);          //CD_epsilon
//           fprintf(outf, "\t%.3f", latmin[i][j]);                        //lateral gap min
//           fprintf(outf, "\t%.3f", radmin[i][j]);                        //radial gap min
//           fprintf(outf, "\t%.3f", offmin[i][j]);                        //offset min
//           fprintf(outf, "\t%.3f", L + radmin[i][j] - offmin[i][j]);     //D-periodicity min
//           fprintf(outf, "\t%.3f", eCDmin[i][j]);                        //E_CD min
//           fprintf(outf, "\t%.3f", eLJmin[i][j]);                        //E_LJ min
//           fprintf(outf, "\t%.3f", eTotmin[i][j]);                       //E_tot min
//           if(xyz_outputs){
//             string xyzfile = xyzbasefile;
//             xyzfile += "-"
//               + replace_char(clean_to_string(lj_min + (i) * lj_stepsize),'.','_')
//               + "-"
//               +  replace_char(clean_to_string(cd_min + (j) * cd_stepsize),'.','_')
//               + ".xyz";
//             // cout << xyzfile << endl;
//             output_xyz(xyzfile,3,latmin[i][j],radmin[i][j],offmin[i][j]);
//           }
//       }
//       fprintf(outf, "\n");
//   }
//   fclose(outf);
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
