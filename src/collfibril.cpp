#include "main.hpp"
#include "collfibril.hpp"

/* Variables */
extern parameters_ parameters;
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

void collagenFibril::singleEmin()
{
  double latGap_ = mol.diameterAtom;
  double radGap_ = mol.diameterAtom;
  double offset_ = mol.diameterAtom;
  double energy_ = 1e16;

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
          // if(xyz_outputs){
          //   string xyzfile = xyzbasefile;
          //   xyzfile += "-"
          //     + replace_char(clean_to_string(lj_min + (i) * lj_stepsize),'.','_')
          //     + "-"
          //     +  replace_char(clean_to_string(cd_min + (j) * cd_stepsize),'.','_')
          //     + ".xyz";
          //   // cout << xyzfile << endl;
          //   output_xyz(xyzfile,3,latmin[i][j],radmin[i][j],offmin[i][j]);
          // }
      }
      fprintf(outf, "\n");
  }
  fclose(outf);
}
