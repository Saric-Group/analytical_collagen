#include "main.hpp"
#include "parse.hpp"
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
  std::string xyzbasefile = "";

  if(flags.xyz_outputs){
    xyzbasefile = filePaths.outputpath;
  }

  filePaths.outputpath += ".dat";

  /* Write header for outputfile */
  header(filePaths.outputpath);

  /* Write results to file */
  FILE *outf;
  outf = fopen(filePaths.outputpath.c_str(), "a");
  /* Both LJ and CD potential */
  for (int i = 0; i < parameters.lj_steps; i++) {
      for (int j = 0; j < parameters.cd_steps; j++) {
          fprintf(outf, "\n");
          fprintf(outf, "%.3f", parameters.lj_min + (i) * parameters.lj_stepsize);            //LJ_epsilon
          fprintf(outf, "\t%.3f", parameters.cd_min + (j) * parameters.cd_stepsize);          //CD_epsilon
          fprintf(outf, "\t%.3f", latmin[i][j]);                        //lateral gap min
          fprintf(outf, "\t%.3f", radmin[i][j]);                        //radial gap min
          fprintf(outf, "\t%.3f", offmin[i][j]);                        //offset min
          fprintf(outf, "\t%.3f", mol.length + radmin[i][j] - offmin[i][j]);     //D-periodicity min
          fprintf(outf, "\t%.3f", eCDmin[i][j]);                        //E_CD min
          fprintf(outf, "\t%.3f", eLJmin[i][j]);                        //E_LJ min
          fprintf(outf, "\t%.3f", eTotmin[i][j]);                       //E_tot min
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

  /* We do not vary the lateral gap here */
  for (int rad = 0; rad <= mol.numAtoms; rad++) {
    // std::cout << "\nrad: " << rad;
    offset_ = radGap_;
    for (int off = rad; off <= mol.numAtoms; off++) {
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
      offset_ += mol.distanceAtoms;
    }
    radGap_ += mol.distanceAtoms;
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
