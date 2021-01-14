#include <chrono>           /* Time */
#include <iostream>         /* cout, ... */
#include <fstream>          /* ifstream, ... */
#include <math.h>           /* ceil, ... */
#include <vector>           /* vector, ... */

using namespace std;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;

/* Variables */
string file;

int layers = 2;             /* Number of layers including focus layer */

double distance_atoms;      /* distance between neighbouring atoms */
double diameter_atom = 1.12;    /* diameter of atom */
int N = 0;                  /* number of atoms in one molecule */
double L;                   /* length of collagen molecule */
double box;                 /* length of periodic box */
double cd_cutoff = 2.0;     /* cutoff of cd potential */
double lj_cutoff = 2.0;     /* cutoff of lj potential */
double max_cutoff = max(cd_cutoff, lj_cutoff);


  // argument overrides
  // takes -i and -o as inputs

  // -i input file
  // -o output file

string inputpath = "./charge_distribution";
string outputpath = "./energy_min.dat";
string file_extension = ".dat";
double dist_atoms = 0.255;

bool set_output = false;

bool charge_hashed_outputs = false;



vector<double> charge, charge_ind;  /* vector for the charges of each atom */
vector<int> type, type_ind; /* vector for the type of each atom */

/* Functions */
void readAtoms(string &file);
void readTypes(string &file);

double distance(double pos, double first, int n, double lat_gap);

double factorLJ(double d);
double LJ_per_mol(double pos, double dx, double ref, double lat_gap);
double newLJ_per_mol(double pos, double dx, double ref, double lat_gap);
double LJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap);
double newLJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap);
double compLJ(double lat_gap, double rad_gap, double offset, int layers);
double newLJ(double lat_gap, double rad_gap, double offset, int layers);

double factorCD(double q1, double q2, double d);
double CD_per_mol(double pos, double q1, double dx, double ref,
                  double lat_gap);
double CD_layer(double pos, double q1, double dx, int layer, double offset,
                double lat_gap);
double compCD(double lat_gap, double rad_gap, double offset, int layers);

void header(string &file);
void singleEmin(string &file, int layers);
void multipleEmin(string &file, int layers, int numbers);
/* map functions */
double energy(double lat_gap, double rad_gap, double offset, int layer1,
              int layer2, double lj, double cd, bool lj_on);
double intra_layer_energy(double lat_gap, double rad_gap, double offset,
                          double lj, double cd, bool lj_on);


int main(int argc, char const *argv[])
{
  bool input_override = false;
  bool output_override = false;
  bool distance_overrride = false;

  if(argc > 0)
  {
    for(int i=0; i<argc; ++i)
      {
        string strarg = argv[i];
        if(input_override)
        {
          inputpath = argv[i];
          input_override = false;
        }
        if(output_override)
        {
          outputpath = argv[i];
          output_override = false;
        }
        if(distance_overrride)
        {
          dist_atoms = stod(argv[i]);
          distance_overrride = false;
        }
        if(strarg.compare("-i")==0)
        {
          input_override = true;
        }
        if(strarg.compare("-o")==0)
        {
          output_override = true;
        }
        if(strarg.compare("-d")==0)
        {
          distance_overrride = true;
        }
        if(strarg.compare("-c")==0)
        {
          charge_hashed_outputs = true;
        }
      }
  }


  time_point start = chrono::high_resolution_clock::now();
  cout << ":: Computation start";

  /************************************************************************/

  file = inputpath;
  // file = "./dis/charge_distribution_36";            /* Anne's collagen */
  //file = "./dis/charge_distribution_1054";       /* Native collagen */
  distance_atoms = dist_atoms;
  // distance_atoms = 0.255;   /* Anne's length */
  //distance_atoms = 0.285;   /* to account for length of 300nm */
  readAtoms(file);
  string file2 = "./dis/atom_types_1054";
  readTypes(file2);
  L = (N - 1) * distance_atoms;
  /************************************************************************/

  /************************************************************************/
  cout << "\n\n:: running computation...";
  /* Only global Emin */

  file=outputpath;

  singleEmin(file, layers);

  /* First numbers Emin */
  // int numbers = 10;
  // file = "./data/multiple/";
  // multipleEmin(file, layers, numbers);

  /* Individual part for maps */
  // bool lj_on = false;
  // double rad_gap, offset, cell_energy, cell_energy_min, min_r, min_o;
  // cell_energy_min = 1e10;
  // for (int rad = 1; rad <= N; rad++) {
  //     cout << "\n --> " << (100. * rad) / (1. * N);
  //     cout << "\% done...";
  //   rad_gap = rad * distance_atoms;
  //   for (int off = 0; off <= N; off++) {
  //     offset = (rad + off) * distance_atoms;
  //     /* Here calculation of energies for the above given geometries */
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 0, 1, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 0, 2, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 1, 2, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 1, 3, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 2, 3, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 2, 4, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 3, 4, 1.0, 14.0, lj_on);
  //     cell_energy += energy(diameter_atom, rad_gap, offset, 3, 5, 1.0, 14.0, lj_on);
  //
  //     cell_energy += 5.0 * intra_layer_energy(diameter_atom, rad_gap, offset, 1.0, 14.0, lj_on);
  //
  //     if (cell_energy <= cell_energy_min) {
  //       cell_energy_min = cell_energy;
  //       min_r = rad_gap;
  //       min_o = offset;
  //     }
  //   }
  // }
  //
  // cout << "\n\nResult:";
  // cout << "\ncell_enery = " << cell_energy_min;
  // cout << "\nmin_r = " << min_r;
  // cout << "\nmin_o = " << min_o;

  // rad_gap = 34.0;
  // offset = 67.0;
  // cout << "\nenergy = " << energy(diameter_atom, rad_gap, offset, 0, 2, 1.0, 14.0, lj_on);
  // cout << "\nintra_energy = " << 5.0 * intra_layer_energy(diameter_atom, rad_gap, offset, 1.0, 14.0, lj_on);

  /* Individual part for maps */


  /* Individual part for energy ratio */
  // file = "./data/hydro/paramspace.dat";
  // FILE *outf;
  // outf = fopen(file.c_str(), "w");
  // fprintf(outf, "#atoms per molecule N = %i", N);
  // fprintf(outf, "\n#distance_atoms = %.3f", distance_atoms);
  // fprintf(outf, "\n#molecule_length = %.3f", L);
  // fprintf(outf, "\n#layers = %i", layers);
  // // fprintf(outf, "\n\n#ratio_CD_LJ");
  // fprintf(outf, "\n\n#CD_param");
  // fprintf(outf, "\tradial_gap_min");
  // fprintf(outf, "\toffset_min");
  // fprintf(outf, "\tE_CD_min");
  // fprintf(outf, "\tE_LJ_min");
  // fprintf(outf, "\tE_total_min");
  // double cd, lj, rad_gap, offset, help;
  // vector<double> radmin, offmin, eCDmin, eLJmin, eTotmin;
  // radmin.assign(1000, 0);
  // offmin.assign(1000, 0);
  // eLJmin.assign(1000, 1e6);
  // eCDmin.assign(1000, 1e6);
  // eTotmin.assign(1000, 1e6);
  // for (int rad = 1; rad <= N; rad++) {
  //   rad_gap = rad * distance_atoms;
  //   cout << "\n --> " << (100. * rad) / (1. * N);
  //   cout << "\% done...";
  //   for (int off = 0; off <= N - rad; off++) {
  //     offset = (rad + off) * distance_atoms;
  //     lj = newLJ(diameter_atom, rad_gap, offset, layers);
  //     // lj = 0;
  //     cd = compCD(diameter_atom, rad_gap, offset, layers);
  //     for (int i = 0; i < 1000; i++) {
  //       help = cd * (0.001 + i * 0.001);
  //       if (lj + help < eTotmin[i]) {
  //         eLJmin[i] = lj;
  //         eCDmin[i] = help;
  //         eTotmin[i] = lj + help;
  //         radmin[i] = rad_gap;
  //         offmin[i] = offset;
  //       }
  //     }
  //   }
  // }
  // for (int i = 0; i < 1000; i++) {
  //   fprintf(outf, "\n");
  //   fprintf(outf, "%.3f", (0.001 + i * 0.001));
  //   fprintf(outf, "\t%.3f", radmin[i]);
  //   fprintf(outf, "\t%.3f", offmin[i]);
  //   fprintf(outf, "\t%.3f", eCDmin[i]);
  //   fprintf(outf, "\t%.3f", eLJmin[i]);
  //   fprintf(outf, "\t%.3f", eTotmin[i]);
  // }
  // fclose(outf);
  // cout << "\nGlobal minimum energy = " << mini;
  // cout << " at r = " << mini_r << " and o = " << mini_o;
  /* Individual part for energy ratio */
  /************************************************************************/

  /************************************************************************/
  time_point end = chrono::high_resolution_clock::now();
  chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
  cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";
  /************************************************************************/

  return 0;
}


/* Functions */
void readAtoms(string &file)
{
    cout << "\n\n:: reading atom configuration...";
    string line;
    ifstream myfile;
    myfile.open(file.c_str(), ios::in);
    while (myfile.peek() != EOF) {
        getline(myfile, line);
        N++;
        charge.push_back(atof(line.c_str()));
        if (charge[N - 1] != 0) {
          charge_ind.push_back(N - 1);
        }
        // cout << "\natom " << N << " has charge " << line;
    }
    myfile.close();
    cout << "\n -> " << N << " atoms read.";
    L = (N - 1) * distance_atoms;
    cout << " Molecule length L = " << L << ".";
}

void readTypes(string &file)
{
  int counter = 0;
  string line;
  ifstream myfile;
  myfile.open(file.c_str(), ios::in);
  while (myfile.peek() != EOF) {
    getline(myfile, line);
    if (line.at(0) != '#') {
      counter++;
      type.push_back(atoi(line.c_str()));
      if (type[counter - 1] == 2 ||
          type[counter - 1] == 10 ||
          type[counter - 1] == 13 ||
          type[counter - 1] == 15 ||
          type[counter - 1] == 20) {
          type_ind.push_back(counter - 1);
      }
    }
  }
  cout << "\n -> " << type_ind.size() << " hydrophobic atoms found.";
}

double distance(double pos, double first, int n, double lat_gap)
{
    double d;
    d = first + n * distance_atoms - pos;
    d *= d;
    d += lat_gap * lat_gap;
    d = sqrt(d);
    return d;
}

double factorLJ(double d)
{
    return 4. * (pow(1. / d, 12.) - pow(1. / d, 6.));
}

double LJ_per_mol(double pos, double dx, double ref, double lat_gap)
{
  double d, sum;
  int left, right;
  bool first = true;
  sum = 0;

  left = ceil(((pos - dx) - ref) / distance_atoms);
  right = floor(((pos + dx) - ref) / distance_atoms);

  if (left < N && right >= 0) {
      /* left end */
      if (left < 0 && right < N) {
          /* left is replaced by 0 here */
          if (first) {
              sum = 0;
              for (int i = 0; i <= right; i++) {
                  d = distance(pos, ref, i, lat_gap);
                  sum += factorLJ(d);
              }
              first = false;
          } else {
              d = distance(pos, ref, 0, lat_gap);
              sum += factorLJ(d);
          }
          if (left == -1) first = true;
      }

      /* bulk */
      if (left >= 0 && right < N) {
          if (first) {
              sum = 0;
              for (int i = 0; i <= right - left; i++) {
                  d = distance(pos, ref, left + i, lat_gap);
                  sum += factorLJ(d);
              }
              first = false;
          }
          if (right == N - 1) first = true;
      }

      /* right end */
      if (left >= 0 && right >= N) {
          /* right is replaced by N - 1 */
          if (first) {
              sum = 0;
              for (int i = 0; i <= N - 1 - left; i++) {
                  d = distance(pos, ref, left + i, lat_gap);
                  sum += factorLJ(d);
              }
              first = false;
          } else {
              d = distance(pos, ref, N, lat_gap);
              sum -= factorLJ(d);
          }
          if (left == N - 1) first = true;
      }
  }
  return sum;
}

double newLJ_per_mol(double pos, double dx, double ref, double lat_gap)
{
  double d, sum;
  int left, right, atomtype;
  sum = 0;

  left = ceil(((pos - dx) - ref) / distance_atoms);
  right = floor(((pos + dx) - ref) / distance_atoms);

  if (left < N && right >= 0) {
      left = max(left, 0);
      right = min(right, N - 1);
      for (int i = 0; i <= right - left; i++) {
          atomtype = type[left + i];
          if (atomtype == 2 ||
              atomtype == 10 ||
              atomtype == 13 ||
              atomtype == 15 ||
              atomtype == 20) {
              d = distance(pos, ref, left + i, lat_gap);
              // cout << "\nhit with " << factorCD(q1, q2, d) / 150;
              // cout << " to change sum = " << sum / 150;
              sum += factorLJ(d);
              // cout << " to sum = " << sum / 150;
          }
      }
  }
  return sum;
}

double LJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap)
{
  double sum = 0;
  sum += LJ_per_mol(pos, dx, layer * offset - box, lat_gap);
  sum += LJ_per_mol(pos, dx, layer * offset, lat_gap);
  sum += LJ_per_mol(pos, dx, layer * offset + box, lat_gap);

  return sum;
}

double newLJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap)
{
  double sum = 0;
  sum += newLJ_per_mol(pos, dx, layer * offset - box, lat_gap);
  sum += newLJ_per_mol(pos, dx, layer * offset, lat_gap);
  sum += newLJ_per_mol(pos, dx, layer * offset + box, lat_gap);

  return sum;
}

double compLJ(double lat_gap, double rad_gap, double offset, int layers)
{
  double pos, sum;
  int left, right;
  double dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
  sum = 0;
  box = L + rad_gap;
  for (int atom = 1; atom <= N; atom++) {
      pos = (atom - 1) * distance_atoms;

      /* Atom atom feels from molecules from layer i: */
      for (int layer = 1; layer < layers; layer++) {
        sum += LJ_layer(pos, dx, layer, offset, lat_gap);
      }


      /* Atom atom feels from left molecule: */
      left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
      right = floor(((pos + lj_cutoff) + box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos + box - (left + i) * distance_atoms);
              sum += factorLJ(d);
          }
      }

      /* Atom atom feels from right molecule: */
      left = ceil(((pos - lj_cutoff) - box) / distance_atoms);
      right = floor(((pos + lj_cutoff) - box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos - box - (left + i) * distance_atoms);
              sum += factorLJ(d);
          }
      }
  }
  return sum;
}

double newLJ(double lat_gap, double rad_gap, double offset, int layers)
{
  int atomtype, left, right;
  double pos, sum, d;
  double dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
  sum = 0;
  box = L + rad_gap;
  for (int atom = 0; atom < (int) type_ind.size(); atom++) {
    pos = type_ind[atom] * distance_atoms;

    /* Atom atom feels from molecules from layer i: */
    for (int layer = 1; layer < layers; layer++) {
      sum += newLJ_layer(pos, dx, layer, offset, lat_gap);
    }

    /* Atom atom feels from left molecule: */
    left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
    right = floor(((pos + lj_cutoff) + box) / distance_atoms);
    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            atomtype = type[left + i];
            if (atomtype == 2 ||
                atomtype == 10 ||
                atomtype == 13 ||
                atomtype == 15 ||
                atomtype == 20) {
                d = abs(pos + box - (left + i) * distance_atoms);
                sum += factorLJ(d);
            }
        }
    }

    /* Atom atom feels from right molecule: */
    left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
    right = floor(((pos + cd_cutoff) - box) / distance_atoms);
    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            atomtype = type[left + i];
            if (atomtype == 2 ||
                atomtype == 10 ||
                atomtype == 13 ||
                atomtype == 15 ||
                atomtype == 20) {
                d = abs(pos - box - (left + i) * distance_atoms);
                sum += factorLJ(d);
            }
        }
    }
  }
  return sum;
}

double factorCD(double q1, double q2, double d)
{
    double tmp = 0;
    double CDconst = 1.0;
    double CDdamping = 1;

    tmp = (CDconst * q1 * q2) / d;
    tmp *= exp(-1.0 * CDdamping * d);

    return tmp;
}

double CD_per_mol(double pos, double q1, double dx, double ref, double lat_gap)
{
    double d, sum, q2;
    int left, right;
    sum = 0;

    left = ceil(((pos - dx) - ref) / distance_atoms);
    right = floor(((pos + dx) - ref) / distance_atoms);

    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            q2 = charge[left + i];
            if (q2 != 0) {
                d = distance(pos, ref, left + i, lat_gap);
                // cout << "\nhit with " << factorCD(q1, q2, d) / 150;
                // cout << " to change sum = " << sum / 150;
                sum += factorCD(q1, q2, d);
                // cout << " to sum = " << sum / 150;
            }
        }
    }
    return sum;
}
double CD_layer(double pos, double q1, double dx, int layer, double offset,
                double lat_gap)
{
  double sum = 0;
  sum += CD_per_mol(pos, q1, dx, layer * offset - box, lat_gap);
  sum += CD_per_mol(pos, q1, dx, layer * offset, lat_gap);
  sum += CD_per_mol(pos, q1, dx, layer * offset + box, lat_gap);

  return sum;
}
double compCD(double lat_gap, double rad_gap, double offset, int layers)
{
  double pos, sum, q1, q2, d;
  int left, right;
  double dx = sqrt(cd_cutoff * cd_cutoff - lat_gap * lat_gap);
  sum = 0;
  box = L + rad_gap;
  for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
      q1 = charge[charge_ind[atom]];
      pos = charge_ind[atom] * distance_atoms;

      /* Atom atom feels from molecules from layer i: */
      for (int layer = 1; layer < layers; layer++) {
        sum += CD_layer(pos, q1, dx, layer, offset, lat_gap);
      }

      /* Atom atom feels from left molecule: */
      left = ceil(((pos - cd_cutoff) + box) / distance_atoms);
      right = floor(((pos + cd_cutoff) + box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              q2 = charge[left + i];
              if (q2 != 0) {
                  d = abs(pos + box - (left + i) * distance_atoms);
                  sum += factorCD(q1, q2, d);
              }
          }
      }

      /* Atom atom feels from right molecule: */
      left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
      right = floor(((pos + cd_cutoff) - box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              q2 = charge[left + i];
              if (q2 != 0) {
                  d = abs(pos - box - (left + i) * distance_atoms);
                  sum += factorCD(q1, q2, d);
              }
          }
      }
  }
  return sum;
}

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

void singleEmin(string &file, int layers)
{
  double lat_gap = diameter_atom;
  double rad_gap, offset;

  // int rad_max = (int) ceil(L / (layers - 1.0) / distance_atoms);
  int rad_max = (int) ceil((L + 2. * max_cutoff) / 0.1);
  int off_max;

  double lj, lje;
  double cd, cde;
  vector<vector<double> > latmin, radmin, offmin, eCDmin, eLJmin, eTotmin;



  /* Assign needed data vectors */
  latmin.assign(50, vector<double> (20, 0));
  radmin.assign(50, vector<double> (20, 0));
  offmin.assign(50, vector<double> (20, 0));
  eLJmin.assign(50, vector<double> (20, 1e6));
  eCDmin.assign(50, vector<double> (20, 1e6));
  eTotmin.assign(50, vector<double> (20, 1e6));

  int lat_max = (int) ceil(max_cutoff / (0.1 * diameter_atom));
  for (int lat = 0; lat <= lat_max; lat++) {
    lat_gap = lat * 0.1 * diameter_atom;
  for (int rad = 0; rad <= rad_max; rad++) {
    // cout << "\n --> " << (100. * rad) / (1. * rad_max);
    // cout << "\% done ...";

    // rad_gap = rad * distance_atoms;
    rad_gap = rad * 0.1;

    // double tmp = (L - (layers - 1.0) * rad_gap) / layers / distance_atoms;
    // double tmp = (L - rad_gap) / distance_atoms;
    // off_max = (int) ceil(tmp);
    off_max = (int) ceil(L / 0.1);

    for (int off = 0; off <= off_max; off++) {
      // offset = rad_gap + off * distance_atoms;
      offset = off * 0.1;
      /* Energy factors calculation for given parameters above */
      lj = compLJ(lat_gap, rad_gap, offset, layers);
      cd = compCD(lat_gap, rad_gap, offset, layers);
      /* Both CD and LJ potential */
      for (int i = 0; i < 50; i++) {
          lje = (i + 1) * 0.01 * lj;
          for (int j = 0; j < 20; j++) {
              cde = cd / ((j + 1) * 10);
              if (lje + cde < eTotmin[i][j]) {
                  eLJmin[i][j] = lje;
                  eCDmin[i][j] = cde;
                  eTotmin[i][j] = lje + cde;
                  latmin[i][j] = lat_gap;
                  radmin[i][j] = rad_gap;
                  offmin[i][j] = offset;
              }
          }
      }
    }
  }
  }


  if(charge_hashed_outputs)
  {
    //create a unique code for the output file

    unsigned long f_charges = 0;

    string fflags = "";

    unsigned long lulen = sizeof(unsigned long)*8;

    int fcount = 0;

    double chargesum=0.0;
    int chargecount = 0;

    for(int c=0;c<charge.size();c++)
    {
      if(charge[c]!=0){
        chargesum+=charge[c];
        chargecount++;
      }
    }

    double chargemean = chargesum/((double) chargecount);

    fflags+="-";
    fflags+=to_string(N);
    fflags+="-";
    fflags+=to_string(dist_atoms);
    fflags+="-";
    fflags+=to_string(chargesum);
    fflags+="-";
    fflags+=to_string(chargemean);
    fflags+="-";



    for(int i=0; i<N; i++)
    {
      if(fcount>=lulen){
        fflags += to_string(f_charges);
        fflags+="-";
        f_charges = 0;
        fcount = 0;
      }
      if((int) charge[i] != 0)
      {
        f_charges+=pow(2,(int)(lulen-fcount-1));
      }

      fcount++;
    }

    fflags += to_string(f_charges);


    if (file.find(file_extension) != string::npos)
    {
      file = file.substr(0, file.find(file_extension));
    }

    file+=fflags;

  }

  if (file.find(file_extension) != string::npos)
  {
      file = file.substr(0, file.find(file_extension));
  }

  file+=".dat";

  /* Write header for outputfile */
  header(file);


  FILE *outf;
  outf = fopen(file.c_str(), "a");
  /* Both LJ and CD potential */
  for (int i = 0; i < 50; i++) {
      for (int j = 0; j < 20; j++) {
          fprintf(outf, "\n");
          fprintf(outf, "%.3f", (i + 1) * 0.01);
          fprintf(outf, "\t%.3f", (j + 1) * 10.);
          // fprintf(outf, "\t%.3f", lat_gap);
          fprintf(outf, "\t%.3f", latmin[i][j]);
          fprintf(outf, "\t%.3f", radmin[i][j]);
          fprintf(outf, "\t%.3f", offmin[i][j]);
          fprintf(outf, "\t%.3f", L + radmin[i][j] - offmin[i][j]);
          fprintf(outf, "\t%.3f", eCDmin[i][j]);
          fprintf(outf, "\t%.3f", eLJmin[i][j]);
          fprintf(outf, "\t%.3f", eTotmin[i][j]);
      }
      fprintf(outf, "\n");
  }
  fclose(outf);
}

/* this function is not maintained anymore */
void multipleEmin(string &file, int layers, int number)
{
  double lat_gap = diameter_atom;
  double rad_gap, offset;

  int rad_max = (int) ceil(L / (layers - 1.0) / distance_atoms);
  int off_max;
  double lj, lje;
  double cd, cde;
  vector<vector<vector<double> > > latmin, radmin, offmin;
  vector<vector<vector<double> > > eCDmin, eLJmin, eTotmin;

  /* Write header for outputfile */
  string tmps = file;
  for (int i = 1; i <= number; i++) {
    file = tmps + "Emin_" + to_string(i) + ".dat";
    header(file);
  }

  /* Assign needed data vectors */
  // latmin.assign(10, vector<vector<double> > (50, vector<double> (20, 0)));
  radmin.assign(number, vector<vector<double> > (50, vector<double> (20, 0)));
  offmin.assign(number, vector<vector<double> > (50, vector<double> (20, 0)));
  eLJmin.assign(number, vector<vector<double> > (50, vector<double> (20, 1e6)));
  eCDmin.assign(number, vector<vector<double> > (50, vector<double> (20, 1e6)));
  eTotmin.assign(number, vector<vector<double> > (50,
                                                  vector<double> (20, 1e6)));

  for (int rad = 0; rad <= rad_max; rad++) {
    cout << "\n --> " << (100. * rad) / (1. * rad_max);
    cout << "\% done ...";
    rad_gap = rad * distance_atoms;
    double tmp = (L - (layers - 1.0) * rad_gap) / layers / distance_atoms;
    off_max = (int) ceil(tmp);
    for (int off = 0; off <= off_max; off++) {
      offset = rad_gap + off * distance_atoms;
      /* Energy factors calculation for given parameters above */
      lj = compLJ(lat_gap, rad_gap, offset, layers);
      cd = compCD(lat_gap, rad_gap, offset, layers);
      /* Both CD and LJ potential */
      for (int i = 0; i < 50; i++) {
          lje = (i + 1) * 0.01 * lj / 5;
          for (int j = 0; j < 20; j++) {
              cde = cd / ((j + 1) * 10 / 2);

              for (int k = 0; k < number; k++) {
                if (lje + cde < eTotmin[k][i][j]) {
                    eLJmin[k][i][j] = lje;
                    eCDmin[k][i][j] = cde;
                    eTotmin[k][i][j] = lje + cde;
                    // latmin[k][i][j] = lat_gap;
                    radmin[k][i][j] = rad_gap;
                    offmin[k][i][j] = offset;
                    break;
                }
              }
          }
      }
    }
  }

  FILE *outf;
  for (int l = 1; l <= number; l++) {
    int k = l - 1;
    file = tmps + "Emin_" + to_string(l) + ".dat";
    outf = fopen(file.c_str(), "a");
    /* Both LJ and CD potential */
    for (int i = 0; i < 50; i++) {
        for (int j = 0; j < 20; j++) {
            fprintf(outf, "\n");
            fprintf(outf, "%.3f", (i + 1) * 0.01 / 5.);
            fprintf(outf, "\t%.3f", (j + 1) * 10. / 2.);
            fprintf(outf, "\t%.3f", lat_gap);
            fprintf(outf, "\t%.3f", radmin[k][i][j]);
            fprintf(outf, "\t%.3f", offmin[k][i][j]);
            fprintf(outf, "\t%.3f", L + radmin[k][i][j] - offmin[k][i][j]);
            fprintf(outf, "\t%.3f", eCDmin[k][i][j]);
            fprintf(outf, "\t%.3f", eLJmin[k][i][j]);
            fprintf(outf, "\t%.3f", eTotmin[k][i][j]);
        }
        fprintf(outf, "\n");
    }
  }
  fclose(outf);
}

double energy(double lat_gap, double rad_gap, double offset, int layer1,
              int layer2, double lj, double cd, bool lj_on)
{
  if (layer1 >= layer2) {
    cout << "error";
    return 0.0;
  }

  double pos, sum, q1;
  double dx = sqrt(cd_cutoff * cd_cutoff - lat_gap * lat_gap);
  sum = 0.0;
  box = L + rad_gap;
  /* cycle through atoms of layer1 */
  for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
    q1 = charge[charge_ind[atom]];
    pos = charge_ind[atom] * distance_atoms;
    /* focus atom1 feels from layer2 */
    sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset - box, lat_gap);
    sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset, lat_gap);
    sum += CD_per_mol(pos, q1, dx, (layer2 - layer1) * offset + box, lat_gap);
  }
  sum *= cd;

  /* turn on/off for hydrophobic interactions */
  if (lj_on) {
    dx = sqrt(lj_cutoff * lj_cutoff - lat_gap * lat_gap);
    for (int atom = 1; atom <= N; atom++) {
      pos = (atom - 1) * distance_atoms;
      sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset - box, lat_gap);
      sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset, lat_gap);
      sum += lj * LJ_per_mol(pos, dx, (layer2 - layer1) * offset + box, lat_gap);
    }
  }
  return sum;
}

double intra_layer_energy(double lat_gap, double rad_gap, double offset,
                          double lj, double cd, bool lj_on)
{
  double pos, sum, q1, q2, d;
  int left, right;
  sum = 0;
  box = L + rad_gap;
  for (int atom = 0; atom < (int) charge_ind.size(); atom++) {
    q1 = charge[charge_ind[atom]];
    pos = charge_ind[atom] * distance_atoms;

    /* Atom atom feels from left molecule: */
    left = ceil(((pos - cd_cutoff) + box) / distance_atoms);
    right = floor(((pos + cd_cutoff) + box) / distance_atoms);
    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            q2 = charge[left + i];
            if (q2 != 0) {
                d = abs(pos + box - (left + i) * distance_atoms);
                sum += factorCD(q1, q2, d);
            }
        }
    }

    /* Atom atom feels from right molecule: */
    left = ceil(((pos - cd_cutoff) - box) / distance_atoms);
    right = floor(((pos + cd_cutoff) - box) / distance_atoms);
    if (left < N && right >= 0) {
        left = max(left, 0);
        right = min(right, N - 1);
        for (int i = 0; i <= right - left; i++) {
            q2 = charge[left + i];
            if (q2 != 0) {
                d = abs(pos - box - (left + i) * distance_atoms);
                sum += factorCD(q1, q2, d);
            }
        }
    }
  }
  sum *= cd;

  if (lj_on) {
    for (int atom = 1; atom <= N; atom++) {
      /* Atom atom feels from left molecule: */
      left = ceil(((pos - lj_cutoff) + box) / distance_atoms);
      right = floor(((pos + lj_cutoff) + box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos + box - (left + i) * distance_atoms);
              sum += lj * factorLJ(d);
          }
      }

      /* Atom atom feels from right molecule: */
      left = ceil(((pos - lj_cutoff) - box) / distance_atoms);
      right = floor(((pos + lj_cutoff) - box) / distance_atoms);
      if (left < N && right >= 0) {
          left = max(left, 0);
          right = min(right, N - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos - box - (left + i) * distance_atoms);
              sum += lj * factorLJ(d);
          }
      }
    }
  }

  return sum;
}
