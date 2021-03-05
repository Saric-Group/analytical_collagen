#include "collfibril.hpp"

/* Variables */
extern filePaths_ filePaths;
extern flags_ flags;

/* Functions */
double collagenFibril::distance(double pos, double first, int n, double lat_gap)
{
    double d;
    d = first + n * mol.distanceAtoms - pos;
    d *= d;
    d += lat_gap * lat_gap;
    d = sqrt(d);
    return d;
}

double collagenFibril::factorLJ(double d)
{
    return 4. * (pow(1. / d, 12.) - pow(1. / d, 6.));
}

double collagenFibril::LJ_per_mol(double pos, double dx, double ref, double lat_gap)
{
  double d, sum;
  int left, right;
  bool first = true;
  sum = 0;

  left = ceil(((pos - dx) - ref) / mol.distanceAtoms);
  right = floor(((pos + dx) - ref) / mol.distanceAtoms);

  if (left < mol.numAtoms && right >= 0) {
      /* left end */
      if (left < 0 && right < mol.numAtoms) {
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
      if (left >= 0 && right < mol.numAtoms) {
          if (first) {
              sum = 0;
              for (int i = 0; i <= right - left; i++) {
                  d = distance(pos, ref, left + i, lat_gap);
                  sum += factorLJ(d);
              }
              first = false;
          }
          if (right == mol.numAtoms - 1) first = true;
      }

      /* right end */
      if (left >= 0 && right >= mol.numAtoms) {
          /* right is replaced by N - 1 */
          if (first) {
              sum = 0;
              for (int i = 0; i <= mol.numAtoms - 1 - left; i++) {
                  d = distance(pos, ref, left + i, lat_gap);
                  sum += factorLJ(d);
              }
              first = false;
          } else {
              d = distance(pos, ref, mol.numAtoms, lat_gap);
              sum -= factorLJ(d);
          }
          if (left == mol.numAtoms - 1) first = true;
      }
  }
  return sum;
}

double collagenFibril::LJ_layer(double pos, double dx, int layer, double offset,
                double lat_gap, double box)
{
  double sum = 0;
  sum += LJ_per_mol(pos, dx, layer * offset - box, lat_gap);
  sum += LJ_per_mol(pos, dx, layer * offset, lat_gap);
  sum += LJ_per_mol(pos, dx, layer * offset + box, lat_gap);

  return sum;
}

double collagenFibril::compLJ(double lat_gap, double rad_gap, double offset)
{
  double pos, sum, box;
  int left, right;
  double dx = sqrt(parameters.lj_cutoff * parameters.lj_cutoff - lat_gap * lat_gap);
  sum = 0;
  box = mol.length + rad_gap;
  for (int atom = 1; atom <= mol.numAtoms; atom++) {
      pos = (atom - 1) * mol.distanceAtoms;

      /* Atom atom feels from molecules from layer i: */
      for (int layer = 1; layer < layers; layer++) {
        sum += LJ_layer(pos, dx, layer, offset, lat_gap, box);
      }


      /* Atom atom feels from left molecule: */
      left = ceil(((pos - parameters.lj_cutoff) + box) / mol.distanceAtoms);
      right = floor(((pos + parameters.lj_cutoff) + box) / mol.distanceAtoms);
      if (left < mol.numAtoms && right >= 0) {
          left = std::max(left, 0);
          right = std::min(right, mol.numAtoms - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos + box - (left + i) * mol.distanceAtoms);
              sum += factorLJ(d);
          }
      }

      /* Atom atom feels from right molecule: */
      left = ceil(((pos - parameters.lj_cutoff) - box) / mol.distanceAtoms);
      right = floor(((pos + parameters.lj_cutoff) - box) / mol.distanceAtoms);
      if (left < mol.numAtoms && right >= 0) {
          left = std::max(left, 0);
          right = std::min(right, mol.numAtoms - 1);
          for (int i = 0; i <= right - left; i++) {
              double d = abs(pos - box - (left + i) * mol.distanceAtoms);
              sum += factorLJ(d);
          }
      }
  }
  return sum;
}

double collagenFibril::factorCD(double q1, double q2, double d)
{
    double tmp = 0;
    double CDconst = 1.0;
    double CDdamping = 1;

    tmp = (CDconst * q1 * q2) / d;
    tmp *= exp(-1.0 * CDdamping * d);

    return tmp;
}

double collagenFibril::CD_per_mol(double pos, double q1, double dx, double ref, double lat_gap)
{
    double d, sum, q2;
    int left, right;
    sum = 0;

    left = ceil(((pos - dx) - ref) / mol.distanceAtoms);
    right = floor(((pos + dx) - ref) / mol.distanceAtoms);

    if (left < mol.numAtoms && right >= 0) {
        left = std::max(left, 0);
        right = std::min(right, mol.numAtoms - 1);
        for (int i = 0; i <= right - left; i++) {
            q2 = mol.charges[left + i];
            if (q2 != 0) {
                d = distance(pos, ref, left + i, lat_gap);
                sum += factorCD(q1, q2, d);
            }
        }
    }
    return sum;
}

double collagenFibril::CD_layer(double pos, double q1, double dx, int layer, double offset,
                double lat_gap, double box)
{
  double sum = 0;
  sum += CD_per_mol(pos, q1, dx, layer * offset - box, lat_gap);
  sum += CD_per_mol(pos, q1, dx, layer * offset, lat_gap);
  sum += CD_per_mol(pos, q1, dx, layer * offset + box, lat_gap);

  return sum;
}

double collagenFibril::compCD(double lat_gap, double rad_gap, double offset)
{
  double pos, sum, q1, q2, d, box;
  int left, right;
  double dx = sqrt(parameters.cd_cutoff * parameters.cd_cutoff - lat_gap * lat_gap);
  sum = 0;
  box = mol.length + rad_gap;
  for (int atom = 0; atom < (int) mol.nonZeroChargesIndex.size(); atom++) {
      q1 = mol.charges[mol.nonZeroChargesIndex[atom]];
      pos = mol.nonZeroChargesIndex[atom] * mol.distanceAtoms;

      /* Atom atom feels from molecules from layer i: */
      for (int layer = 1; layer < layers; layer++) {
        sum += CD_layer(pos, q1, dx, layer, offset, lat_gap, box);
      }

      /* Atom atom feels from left molecule: */
      left = ceil(((pos - parameters.cd_cutoff) + box) / mol.distanceAtoms);
      right = floor(((pos + parameters.cd_cutoff) + box) / mol.distanceAtoms);
      if (left < mol.numAtoms && right >= 0) {
          left = std::max(left, 0);
          right = std::min(right, mol.numAtoms - 1);
          for (int i = 0; i <= right - left; i++) {
              q2 = mol.charges[left + i];
              if (q2 != 0) {
                  d = abs(pos + box - (left + i) * mol.distanceAtoms);
                  sum += factorCD(q1, q2, d);
              }
          }
      }

      /* Atom atom feels from right molecule: */
      left = ceil(((pos - parameters.cd_cutoff) - box) / mol.distanceAtoms);
      right = floor(((pos + parameters.cd_cutoff) - box) / mol.distanceAtoms);
      if (left < mol.numAtoms && right >= 0) {
          left = std::max(left, 0);
          right = std::min(right, mol.numAtoms - 1);
          for (int i = 0; i <= right - left; i++) {
              q2 = mol.charges[left + i];
              if (q2 != 0) {
                  d = abs(pos - box - (left + i) * mol.distanceAtoms);
                  sum += factorCD(q1, q2, d);
              }
          }
      }
  }
  return sum;
}

int collagenFibril::produce_xyz(int numatoms, int rows, double distatom, double latgap, double radgap, double offset, int id[], int type[], double charges[], double xpos[], double ypos[], double zpos[])
{
  int natoms = layers * rows * numatoms;

  //cout << natoms << endl;
  //cout << "Atoms. Timestep: 0" << endl;

  std::map<int,int> typemap;
  int typecount = 0;
  for (int i = 0; i < numatoms; i++){
    if(not typemap.count(mol.charges[i])){
      typemap[mol.charges[i]] = typecount;
      typecount += 1;
    }
  }

  for (int l = 0; l < layers; l++){
    for (int r = 0; r < rows; r++){
      for (int a = 0; a < numatoms; a++){
        int i = l * rows * numatoms + r * numatoms + a;
        double drows = (double) rows;
        double dr =(double) r;
        double da =(double) a;
        double dl =(double) l;
        id[i] = i;
        charges[i] = mol.charges[a];
        type[i] = typemap[mol.charges[a]];
        xpos[i] = (((double)numatoms)*dr*distatom) + (dr*radgap) + dl*offset + da*distatom;
        ypos[i] = -dl*latgap;
        zpos[i] = 0.0;
        //cout << id[i] << "  " << xpos[i] << "  " << ypos[i] << "  " << zpos[i] << "  " << type[i] << "  " << charges[i] << endl;
      }
    }
  }

  return 1;

}

int collagenFibril::output_xyz(std::string filename, int rows, double latgap, double radgap, double offset)
{

  int natoms = layers * rows * mol.numAtoms;
  int * id = new int[natoms];
  int * type = new int[natoms];
  double * charges = new double[natoms];
  double * xpos = new double[natoms];
  double * ypos = new double[natoms];
  double * zpos = new double[natoms];

  int err = produce_xyz(mol.numAtoms , rows, mol.distanceAtoms, latgap, radgap, offset, id, type, charges, xpos, ypos, zpos);

  if(err != 1){
    return 0;
  }

  FILE *xyzf;
  xyzf = fopen(filename.c_str(), "w");

  fprintf(xyzf, "%i\n", natoms);
  fprintf(xyzf, "Atoms. Timestep: 0\n");
  for(int i=0; i<natoms;i++){
    std::string xyzrow =
    std::to_string(id[i])
    + "  "
    + clean_to_string(xpos[i])
    + "  "
    + clean_to_string(ypos[i])
    + "  "
    + clean_to_string(zpos[i])
    + "  "
    + std::to_string(type[i])
    + "  "
    + clean_to_string(charges[i])
    + "\n";
    fprintf(xyzf, "%s",xyzrow.c_str());
  }
  fclose (xyzf);

  return err;
}

std::string collagenFibril::hash_output()
{
  //create a unique code for the output file

    unsigned long f_charges = 0;

    std::string fflags = "";

    unsigned long lulen = sizeof(unsigned long)*8;

    int fcount = 0;

    double chargesum=0.0;
    int chargecount = 0;

    for(int c = 0; c < (int) mol.charges.size(); c++)
    {
      if(mol.charges[c] != 0){
        chargesum += mol.charges[c];
        chargecount++;
      }
    }

    double chargemean = chargesum / ((double) chargecount);

    fflags += "-";
    fflags += std::to_string(mol.numAtoms);
    fflags += "-";
    fflags += std::to_string(mol.distanceAtoms);
    fflags += "-";
    fflags += std::to_string(chargesum);
    fflags += "-";
    fflags += std::to_string(chargemean);
    fflags += "-";



    for(int i = 0; i < mol.numAtoms; i++)
    {
      if(fcount >= lulen){
        fflags += std::to_string(f_charges);
        fflags += "-";
        f_charges = 0;
        fcount = 0;
      }
      if((int) mol.charges[i] != 0)
      {
        f_charges += pow(2,(int) (lulen - fcount - 1));
      }

      fcount++;
    }

    fflags += std::to_string(f_charges);


    if (filePaths.outputpath.find(filePaths.file_extension) != std::string::npos)
    {
      filePaths.outputpath = filePaths.outputpath.substr(0, filePaths.outputpath.find(filePaths.file_extension));
    }

    filePaths.outputpath += fflags;
    return filePaths.outputpath;
}

void collagenFibril::header(std::string &file)
{
  FILE *outf;
  outf = fopen(file.c_str(), "w");
  fprintf(outf, "#atoms per molecule N = %i", mol.numAtoms);
  fprintf(outf, "\n#distance_atoms = %.3f", mol.distanceAtoms);
  fprintf(outf, "\n#molecule_length = %.3f", mol.length);
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

void collagenFibril::singleEmin()
{
  std::cout << "\n#\n#";
  double latGap_ = mol.diameterAtom;
  double radGap_ = mol.diameterAtom;
  double offset_ = mol.diameterAtom;

  int latMax = (int) ceil(parameters.max_cutoff / (0.1 * mol.diameterAtom));
  int radMax = (int) ceil((mol.length + 2. * parameters.max_cutoff) / 0.1);
  int offMax;

  double lj, lje;
  double cd, cde;
  std::vector<std::vector<double> > latmin, radmin, offmin, eCDmin, eLJmin, eTotmin;
  /* Assign needed data vectors */
  latmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 0));
  radmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 0));
  offmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 0));
  eLJmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 1e6));
  eCDmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 1e6));
  eTotmin.assign(parameters.lj_steps, std::vector<double> (parameters.cd_steps, 1e6));

  /* Calculate min configurations for different energy settings */
  int total = (latMax + 1) * (radMax + 1);
  int counter = 0;
  for (int lat = 0; lat <= latMax; lat++) {
    latGap_ = lat * 0.1 * mol.diameterAtom;
    for (int rad = 0; rad <= radMax; rad++) {
      radGap_ = rad * 0.1;
      offMax = (int) ceil((mol.length + radGap_) / 0.1);
      for (int off = 0; off <= offMax; off++) {
        offset_ = off * 0.1;
        /* ToDo Energy calculation */
        /***************************/
        lj = compLJ(latGap_, radGap_, offset_);
        cd = compCD(latGap_, radGap_, offset_);
        for (int i = 0; i < parameters.lj_steps; i++) {
            lje = (parameters.lj_min + (i) * parameters.lj_stepsize) * lj;
            for (int j = 0; j < parameters.cd_steps; j++) {
                cde = cd / (parameters.cd_min + (j) * parameters.cd_stepsize);
                if (lje + cde < eTotmin[i][j]) {
                    eLJmin[i][j] = lje;
                    eCDmin[i][j] = cde;
                    eTotmin[i][j] = lje + cde;
                    latmin[i][j] = latGap_;
                    radmin[i][j] = radGap_;
                    offmin[i][j] = offset_;
                }
            }
        }
      }
      counter++;
      // progresssBar(1.0 * counter / total, "Calculating minimized energy configurations");
    }
  }

  if(flags.charge_hashed_outputs)
  {
    hash_output();
  }

  if (filePaths.outputpath.find(filePaths.file_extension) != std::string::npos)
  {
      filePaths.outputpath = filePaths.outputpath.substr(0, filePaths.outputpath.find(filePaths.file_extension));
  }
  if (filePaths.outputpath.find(filePaths.csv_extension) != std::string::npos)
  {
      filePaths.outputpath = filePaths.outputpath.substr(0, filePaths.outputpath.find(filePaths.csv_extension));
  }
  std::string xyzbasefile = "";

  if(flags.xyz_outputs){
    xyzbasefile = filePaths.outputpath;
  }

  FILE *outf;
  if (flags.csv_output) {
    filePaths.outputpath += filePaths.csv_extension;
    outf = fopen(filePaths.outputpath.c_str(), "w");
    fprintf(outf, "LJ_epsilon,");
    fprintf(outf, "CD_epsilon,");
    fprintf(outf, "lateral_gap_min,");
    fprintf(outf, "radial_gap_min,");
    fprintf(outf, "offset_min,");
    fprintf(outf, "D-periodicity_min,");
    fprintf(outf, "E_CD_min,");
    fprintf(outf, "E_LJ_min,");
    fprintf(outf, "E_total_min");
    fclose(outf);
  } else {
    filePaths.outputpath += filePaths.file_extension;
    /* Write header for outputfile */
    header(filePaths.outputpath);
  }

  /* Write results to file */
  outf = fopen(filePaths.outputpath.c_str(), "a");
  /* Both LJ and CD potential */
  std::string separator;
  if (flags.csv_output) {
    separator = ",";
  } else {
    separator = "\t";
  }
  for (int i = 0; i < parameters.lj_steps; i++) {
      for (int j = 0; j < parameters.cd_steps; j++) {
          fprintf(outf, "\n");
          fprintf(outf, "%.3f", parameters.lj_min + (i) * parameters.lj_stepsize);            //LJ_epsilon
          fprintf(outf, "%s%.3f", separator.c_str(), parameters.cd_min + (j) * parameters.cd_stepsize);          //CD_epsilon
          fprintf(outf, "%s%.3f", separator.c_str(), latmin[i][j]);                        //lateral gap min
          fprintf(outf, "%s%.3f", separator.c_str(), radmin[i][j]);                        //radial gap min
          fprintf(outf, "%s%.3f", separator.c_str(), offmin[i][j]);                        //offset min
          fprintf(outf, "%s%.3f", separator.c_str(), mol.length + radmin[i][j] - offmin[i][j]);     //D-periodicity min
          fprintf(outf, "%s%.3f", separator.c_str(), eCDmin[i][j]);                        //E_CD min
          fprintf(outf, "%s%.3f", separator.c_str(), eLJmin[i][j]);                        //E_LJ min
          fprintf(outf, "%s%.3f", separator.c_str(), eTotmin[i][j]);                       //E_tot min
          if(flags.xyz_outputs){
            std::string xyzfile = xyzbasefile;
            xyzfile += "-"
              + replace_char(clean_to_string(parameters.lj_min + (i) * parameters.lj_stepsize),'.','_')
              + "-"
              +  replace_char(clean_to_string(parameters.cd_min + (j) * parameters.cd_stepsize),'.','_')
              + ".xyz";
            // cout << xyzfile << endl;
            output_xyz(xyzfile,3,latmin[i][j],radmin[i][j],offmin[i][j]);
          }
      }
      fprintf(outf, "\n");
  }
  fclose(outf);
}

void collagenFibril::minimizeEnergy()
{
  double latGap_ = mol.diameterAtom;
  double radGap_ = mol.diameterAtom;
  double offset_ = mol.diameterAtom;
  double energy_ = 1e16;
  double a = 1.0;
  int N = 0.5 * (a * mol.numAtoms + 1) * (a * mol.numAtoms + 2);
  int counter = 0;

  /* We do not vary the lateral gap here */
  for (int rad = 0; rad <= a * mol.numAtoms; rad++) {
    // std::cout << "\nrad: " << rad;
    offset_ = radGap_;
    for (int off = rad; off <= a * mol.numAtoms; off++) {
      energy_ = compCD(latGap_, radGap_, offset_);
      // std::cout << "\nradGap: " << radGap_;
      // std::cout << "\toffset: " << offset_;
      // std::cout << "\tenergy: " << energy_;
      if (energy_ < energy) {
        energy = energy_;
        latGap = latGap_;
        radGap = radGap_;
        offset = offset_;
      }
      offset_ += mol.distanceAtoms / a;
      counter++;
      // progresssBar(1.0 * counter / N, " fibril -> minimize energy ");
    }
    radGap_ += mol.distanceAtoms / a;
  }
}

void collagenFibril::writeXYZ()
{
  FILE *outf;
  std::string file = filePaths.outputpath;
  if (file.find(filePaths.file_extension) != std::string::npos)
  {
      file = file.substr(0, file.find(filePaths.file_extension));
  }
  file += ".xyz";
  outf = fopen(file.c_str(), "w");
  int units = 1;
  fprintf(outf, "%i", mol.numAtoms * layers * 3 * units);
  fprintf(outf, "\n");

  double box = mol.length + radGap;
  for (int j = 0; j < units; j++) {
    for (int layer = 0; layer < layers; layer++) {
      /* center molecule */
      for (int i = 0; i < mol.numAtoms; i++) {
        fprintf(outf, "\n%.3f", mol.charges[i]);
        fprintf(outf, "\t%.8f", layer * offset + i * mol.distanceAtoms);
        fprintf(outf, "\t%.8f", (layer + j * layers) * latGap);
        fprintf(outf, "\t%.8f", 0.);
      }
      /* right molecule */
      for (int i = 0; i < mol.numAtoms; i++) {
        fprintf(outf, "\n%.3f", mol.charges[i]);
        fprintf(outf, "\t%.8f", box + layer * offset + i * mol.distanceAtoms);
        fprintf(outf, "\t%.8f", (layer + j * layers) * latGap);
        fprintf(outf, "\t%.8f", 0.);
      }
      /* left molecule */
      for (int i = 0; i < mol.numAtoms; i++) {
        fprintf(outf, "\n%.3f", mol.charges[i]);
        fprintf(outf, "\t%.8f", -box + layer * offset + i * mol.distanceAtoms);
        fprintf(outf, "\t%.8f", (layer + j * layers) * latGap);
        fprintf(outf, "\t%.8f", 0.);
      }
    }
  }
  fclose(outf);
}

void collagenFibril::printMoleculeInfo()
{
  std::cout << "\n#\n#";
  std::cout << "\n# Molecule information:";
  std::cout << "\n#    Number of atoms: " << mol.numAtoms;
  std::cout << "\n#    Number of types: " << mol.numTypes;
  std::cout << "\n#    Molecule length: " << mol.length;
  std::cout << "\n#    Interatomic distance: " << mol.distanceAtoms;
  std::cout << "\n#    Atom diameter: " << mol.diameterAtom;
  std::cout << "\n#    Number of positive charges: " << mol.numPos;
  std::cout << "\n#    Number of negative charges: " << mol.numNeg;
  std::cout << "\n#    Total charge: " << mol.totalCharge;
}

/* Functions  to clean up*/
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
