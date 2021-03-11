#include "md.hpp"


/* Variables */
extern std::random_device device;
extern std::mt19937 generator;
extern filePaths_ filePaths;
extern flags_ flags;


/* Operators */
vec3d operator*(const double &k, const vec3d &a)
{
  return a * k;
}

std::ostream& operator<<(std::ostream &os, const vec3d &a)
{
  os << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return os;
}


/* Functions */
void vec3d::calcLength()
{
  length = *this * *this;
  length = sqrt(length);
}

void vec3d::normalize()
{
  calcLength();
  if (length != 0) {
    *this /= length;
  }
  calcLength();
}

vec3d crossproduct(vec3d a, vec3d b)
{
  vec3d v(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
  return v;
}

vec3d getLinePtA(vec3d centerPoint, double length, double phi, double theta)
{
  double x = centerPoint.x - 0.5 * length * cos(phi) * sin(theta);
  double y = centerPoint.y - 0.5 * length * sin(phi) * sin(theta);
  double z = centerPoint.z - 0.5 * length * cos(theta);
  return vec3d(x, y, z);
}

vec3d getLinePtB(vec3d centerPoint, double length, double phi, double theta)
{
  double x = centerPoint.x + 0.5 * length * cos(phi) * sin(theta);
  double y = centerPoint.y + 0.5 * length * sin(phi) * sin(theta);
  double z = centerPoint.z + 0.5 * length * cos(theta);
  return vec3d(x, y, z);
}

double shortestDistBetweenLines(vec3d a1, vec3d a2, vec3d b1, vec3d b2)
{
  // just Algebra things
  vec3d a = a2 - a1;
  vec3d b = b2 - b1;
  vec3d n = crossproduct(a, b);

  if (n.length == 0) {
    return crossproduct(a, b1 - a1).length / a.length;
  } else {
    double s = b1 * a - a1 * a;
    s *= a * b;
    s /= a * a;
    s -= b1 * b;
    s += a1 * b;
    s /= b * b - (b * a) * (a * b) / (a * a);
    double t = (b1 * a + s * b * a - a1 * a) / (a * a);
    s = std::min(std::max(s, 0.0), 1.0);
    t = std::min(std::max(t, 0.0), 1.0);
    vec3d c1 = a1 + t * a;
    vec3d c2 = b1 + s * b;
    return (c2 - c1).length;
  }
}

bool moleculesOverlap(int a, int b, collagenFibril fib,
                  std::vector<vec3d> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem)
{
  if (b < 0) {  // no molecules exist with index b < 0
    return false;
  }

  double length = fib.mol.length;
  vec3d a_A = getLinePtA(gridPoints[a], length, phi_mem[a], theta_mem[a]);
  vec3d a_B = getLinePtB(gridPoints[a], length, phi_mem[a], theta_mem[a]);
  vec3d b_A = getLinePtA(gridPoints[b], length, phi_mem[b], theta_mem[b]);
  vec3d b_B = getLinePtB(gridPoints[b], length, phi_mem[b], theta_mem[b]);

  double dist = shortestDistBetweenLines(a_A, a_B, b_A, b_B);

  // 100% safety margin -> distance of double atom diameter enforced
  return fib.mol.diameterAtom * 2.0 > dist;
}

bool checkOverlap(int i, int L, std::vector<vec3d> &gridPoints,
                  std::vector<double> &phi_mem,
                  std::vector<double> &theta_mem,
                  collagenFibril fib)
{
  int j;
  /*
   * Max 13 molecules are already placed that need to be considered when
   * placing new molecule i, i.e. the 9 molecules in the layer behind it,
   * three molecules in the column to the left of it,
   * and finally one molecule below it.
   * If only one of those overlaps, we can exit and randomize i's position
   * anew.
   */

  // if second param in moleculesOverlap is < 1 it returns false

  // one molecule below
  if (i % L != 0) { // i % L = 0 are the bottom molecules
    if (moleculesOverlap(i, i - 1, fib, gridPoints, phi_mem, theta_mem)) return true;
  }

  // three molecules to the left
  if (i % (L * L) >= L) { // i % (L * L) < L are the outer left molecules
    j = i - L;
    // molecule center left
    if (moleculesOverlap(i, j, fib, gridPoints, phi_mem, theta_mem)) return true;
    // molecule below left
    if (i % L != 0) {
      if (moleculesOverlap(i, j - 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
    // molecule above left
    if (i % L != L - 1) {
      if (moleculesOverlap(i, j + 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
  }

  // nine molecules in layer behind
  j = i - L * L;
  // molecule center behind
  if (moleculesOverlap(i, j, fib, gridPoints, phi_mem, theta_mem)) return true;
  // molecule below center behind
  if (i % L != 0) {
    if (moleculesOverlap(i, j - 1, fib, gridPoints, phi_mem, theta_mem)) return true;
  }
  // molecule above center behind
  if (i % L != L - 1) {
    if (moleculesOverlap(i, j + 1, fib, gridPoints, phi_mem, theta_mem)) return true;
  }
  if (i % (L * L) >= L) {
    // molecule center left behind
    if (moleculesOverlap(i, j - L, fib, gridPoints, phi_mem, theta_mem)) return true;
    // molecule below left behind
    if (i % L != 0) {
      if (moleculesOverlap(i, j - L - 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
    // molecule above left behind
    if (i % L != L - 1) {
      if (moleculesOverlap(i, j - L + 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
  }
  if (i % (L * L) < L * (L - 1)) { // i % (L * L) >= L * (L - 1) are the outer right molecules
    // molecule center right behind
    if (moleculesOverlap(i, j + L, fib, gridPoints, phi_mem, theta_mem)) return true;
    // molecule below right behind
    if (i % L != 0) {
      if (moleculesOverlap(i, j + L - 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
    // molecule above right behind
    if (i % L != L - 1) {
      if (moleculesOverlap(i, j + L + 1, fib, gridPoints, phi_mem, theta_mem)) return true;
    }
  }

  // placement of molecule i is good, no overlap, return false
  return false;
}

void genTopologyZero(collagenFibril fib, int L)
{
  // number of molecules per sidelength of cubic box
  // #define L 10

  // total number of molecules in box
  int numMol = L * L * L;

  // defines the boundaries of cubic simulation box with length=boxlength
  // to be in line with Anne's volume fraction, we have to calculate:
  double result = 0.0;
  result = (36 - 1) * 0.255 * 1000 / (100 * 100 * 100);
  result /= (fib.mol.numAtoms - 1) * fib.mol.distanceAtoms * numMol;
  result = 1 / result;
  result = pow(result, 1. / 3.);
  // Increased the box ==> lower volume fraction, to keep molecules in box
  // might need to think about boxsize again
	double boxlength = result + 1.0 * fib.mol.length;
	double xlo = -0.5 * boxlength;
	double ylo = -0.5 * boxlength;
	double zlo = -0.5 * boxlength;
	double xhi = 0.5 * boxlength;
	double yhi = 0.5 * boxlength;
	double zhi = 0.5 * boxlength;

  std::string file = filePaths.md_outputpath;
  file += "topology.0time";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  /* LAMMPS Descritption */
  fprintf(outf, "LAMMPS Description\n\n");

  fprintf(outf, "\t%i atoms", numMol * fib.mol.numAtoms);
  fprintf(outf, "\n\t%i bonds", numMol * (fib.mol.numAtoms - 1));
  fprintf(outf, "\n\t%i angles", numMol * (fib.mol.numAtoms - 2));
  fprintf(outf, "\n\t0 dihedrals");
  fprintf(outf, "\n\t0 impropers");

  fprintf(outf, "\n\n\t%i atom types", (int) fib.mol.numTypes);
  fprintf(outf, "\n\t1 bond types");
  fprintf(outf, "\n\t1 angle types");
  fprintf(outf, "\n\t0 dihedral types");
  fprintf(outf, "\n\t0 improper types");

  fprintf(outf, "\n\n%.6f %.6f xlo xhi", xlo, xhi);
  fprintf(outf, "\n%.6f %.6f ylo yhi", ylo, yhi);
  fprintf(outf, "\n%.6f %.6f zlo zhi", zlo, zhi);


  /* Masses */
  fprintf(outf, "\n\nMasses");
  fprintf(outf, "\n");
  for (int i = 0; i < fib.mol.numTypes; i++) {
    fprintf(outf, "\n\t%i 1.000", i + 1);
  }


  /* Atoms */
  cubeGrid grid(boxlength, L);
  std::vector<double> phi_mem, theta_mem;
  fprintf(outf, "\n\nAtoms");
  fprintf(outf, "\n");
  std::uniform_real_distribution<double> dis_real(0.0, M_PI);
  double phi = fib.parametersMD.phi;
  double theta = fib.parametersMD.theta;
  if (!fib.parametersMD.random) {
    phi_mem.assign(numMol, phi);
    theta_mem.assign(numMol, theta);
    if (checkOverlap(numMol / 2, L, grid.gridPoints, phi_mem, theta_mem, fib)) {
      std::cout << "\n#  ERROR: overlapping molecules, choose other orientation angles! ";
      return;
    }
  }
  for (int i = 0; i < numMol; i++) {
    if (fib.parametersMD.random) {
      bool overlap = false;
      do {
        if (overlap) {
          phi_mem.pop_back();
          theta_mem.pop_back();
        }
        phi = 2 * dis_real(generator);
        theta = dis_real(generator);
        phi_mem.push_back(phi);
        theta_mem.push_back(theta);
        overlap = checkOverlap(i, L, grid.gridPoints, phi_mem, theta_mem, fib);
      } while(overlap);
    }
    for (int j = 0; j < fib.mol.numAtoms; j++) {
      fprintf(outf, "\n\t%i", i * fib.mol.numAtoms + j + 1);  // Atom ID
      fprintf(outf, " %i", i + 1);  // Molecule ID
      fprintf(outf, " %i", fib.mol.atomTypes[j]);   // Atom type
      fprintf(outf, " %.6f", fib.mol.charges[j]);   // Atom charge
      fprintf(outf, " %.6f", grid.gridPoints[i].x - (0.5 * fib.mol.length - j * fib.mol.distanceAtoms) * cos(phi) * sin(theta));
      fprintf(outf, " %.6f", grid.gridPoints[i].y - (0.5 * fib.mol.length - j * fib.mol.distanceAtoms) * sin(phi) * sin(theta));
      fprintf(outf, " %.6f", grid.gridPoints[i].z - (0.5 * fib.mol.length - j * fib.mol.distanceAtoms) * cos(theta));
    }
  }

  /* Bonds */
  fprintf(outf, "\n\nBonds");
  fprintf(outf, "\n");
  for (int i = 0; i < numMol; i++) {
    for (int j = 1; j < fib.mol.numAtoms; j++) {
      fprintf(outf, "\n\t%i", i * (fib.mol.numAtoms - 1) + j);  // Bond ID
      fprintf(outf, " 1");    // Bond type
      fprintf(outf, " %i", i * fib.mol.numAtoms + j);   // Atom 1 ID
      fprintf(outf, " %i", i * fib.mol.numAtoms + j + 1);    // Atom 2 ID
    }
  }


  /* Angles */
  fprintf(outf, "\n\nAngles");
  fprintf(outf, "\n");
  for (int i = 0; i < numMol; i++) {
    for (int j = 1; j < fib.mol.numAtoms - 1; j++) {
      fprintf(outf, "\n\t%i", i * (fib.mol.numAtoms - 2) + j);  // Angle ID
      fprintf(outf, " 1");    // Angle type
      fprintf(outf, " %i", i * fib.mol.numAtoms + j);   // Atom 1 ID
      fprintf(outf, " %i", i * fib.mol.numAtoms + j + 1);    // Atom 2 ID
      fprintf(outf, " %i", i * fib.mol.numAtoms + j + 2);    // Atom 3 ID
    }
  }

  fclose(outf);
}

void printScriptVars(FILE *outf, md_var var, int tabs)
{
  fprintf(outf, "\n");
  for (int i = 0; i < tabs; i++) {
    fprintf(outf, "\t");
  }
  fprintf(outf, "for %s in $(seq", var.abr.c_str());

  if (var.size == sizeof(int)) {
    fprintf(outf, " %.i", *(int *) var.vars[0]);
    fprintf(outf, " %.i", *(int *) var.vars[1]);
    fprintf(outf, " %.i)", *(int *) var.vars[2]);
  }
  if (var.size == sizeof(double)) {
    fprintf(outf, " %.*f", var.prec, *(double *) var.vars[0]);
    fprintf(outf, " %.*f", var.prec, *(double *) var.vars[1]);
    fprintf(outf, " %.*f)", var.prec, *(double *) var.vars[2]);
  }

  fprintf(outf, "\n");
  for (int i = 0; i < tabs; i++) {
    fprintf(outf, "\t");
  }
  fprintf(outf, "do");
}

std::vector<md_var> collectMDvars(collagenFibril fib)
{
  std::vector<md_var> md_vars;

  if (fib.parametersMD.kAngle_start + fib.parametersMD.kAngle_inc <= fib.parametersMD.kAngle_end) {
    md_vars.push_back(md_var("ka", sizeof(fib.parametersMD.kAngle_start),
                &fib.parametersMD.kAngle,
                std::vector<double *>{&fib.parametersMD.kAngle_start,
                                      &fib.parametersMD.kAngle_inc,
                                      &fib.parametersMD.kAngle_end}, 1));
  }

  if (fib.parametersMD.LJepsilon_start + fib.parametersMD.LJepsilon_inc <= fib.parametersMD.LJepsilon_end) {
    md_vars.push_back(md_var("lje", sizeof(fib.parametersMD.LJepsilon_start),
                &fib.parametersMD.LJepsilon,
                std::vector<double *>{&fib.parametersMD.LJepsilon_start,
                                      &fib.parametersMD.LJepsilon_inc,
                                      &fib.parametersMD.LJepsilon_end}, 3));
  }

  if (fib.parametersMD.dielectric_start + fib.parametersMD.dielectric_inc <= fib.parametersMD.dielectric_end) {
    md_vars.push_back(md_var("die", sizeof(fib.parametersMD.dielectric_start),
                &fib.parametersMD.dielectric,
                std::vector<double *>{&fib.parametersMD.dielectric_start,
                                      &fib.parametersMD.dielectric_inc,
                                      &fib.parametersMD.dielectric_end}, 1));
  }

  return md_vars;
}

std::string getDir(std::vector<md_var> md_vars, bool pvalues /*= false*/)
{
  std::string dir = "";
  for (int i = 0; i < (int) md_vars.size(); i++) {
    if (i > 0) {
      dir += "_";
    }
    dir += md_vars[i].abr;
    dir += "=";
    if (pvalues) {
      if (md_vars[i].size == sizeof(int)) {
        dir += *(int *) md_vars[i].var;
      } else if (md_vars[i].size == sizeof(double)) {
        std::ostringstream oss;
        oss << std::setprecision(md_vars[i].prec);
        oss << std::fixed;
        oss << *(double *) md_vars[i].var;
        dir += oss.str();
      }
    } else {
      dir += "${" + md_vars[i].abr + "}";
    }
  }
  if (md_vars.size() > 0) {
    dir += "/";
  }
  return dir;
}

std::string getCargs(std::vector<md_var> md_vars)
{
  std::string cargs = "";
  for (int i = 0; i < (int) md_vars.size(); i++) {
    if (i > 0) {
      cargs += " ";
    }
    cargs += "-md_";
    cargs += md_vars[i].abr;
    cargs += " ";
    cargs += "${" + md_vars[i].abr + "}";
  }
  return cargs;
}

void genInSim(collagenFibril fib)
{
  std::string file = filePaths.md_outputpath;
  file += "in.sim";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  std::vector<md_var> md_vars = collectMDvars(fib);
  std::string dir = getDir(md_vars, true);

  fprintf(outf, "log");
  fprintf(outf, "\t./%slog.sim", dir.c_str());

  fprintf(outf, "\n\nunits\tlj");
  fprintf(outf,"\natom_style\tfull");

  fprintf(outf, "\n\nread_data\t");
  fprintf(outf, "./topology.0time");

  fprintf(outf, "\n\nvariable\tLANGEVIN_SEED\tequal\t1337");

  fprintf(outf, "\n\nbond_style\tharmonic");
  fprintf(outf, "\nbond_coeff\t1 500 %.6f", fib.mol.distanceAtoms);

  fprintf(outf, "\n\nangle_style\tharmonic");
  fprintf(outf, "\nangle_coeff\t1 %.1f 180", fib.parametersMD.kAngle);

  fprintf(outf, "\n\npair_style\thybrid/overlay lj/cut %.1f coul/debye 1 %.1f", fib.parametersMD.lj_cutoff, fib.parametersMD.cd_cutoff);
  fprintf(outf, "\ndielectric\t%.1f", fib.parametersMD.dielectric);
  fprintf(outf, "\npair_coeff\t* * lj/cut %.3f 1 %.1f", fib.parametersMD.LJepsilon, fib.parametersMD.lj_cutoff);
  fprintf(outf, "\npair_coeff\t* * coul/debye %.1f", fib.parametersMD.cd_cutoff);
  fprintf(outf, "\npair_modify\tshift yes");
  fprintf(outf, "\nneigh_modify\texclude molecule/intra all");

  fprintf(outf, "\n\nvelocity\tall create 1.0 ${LANGEVIN_SEED}");
  fprintf(outf, "\nvelocity\tall zero angular");
  fprintf(outf, "\nvelocity\tall zero linear");

  fprintf(outf, "\n\ndump\tmydump all custom 10000 ./%s", dir.c_str());
  fprintf(outf, "out_sim.xyz type q x y z");
  fprintf(outf, "\ndump_modify\tmydump sort id");

  if (fib.parametersMD.rigid) {
    fprintf(outf, "\n\nfix\t0 all rigid/nve molecule");
  } else {
    fprintf(outf, "\n\nfix\t0 all nve");
  }
  fprintf(outf, "\nfix\t2 all langevin 1 1 1 ${LANGEVIN_SEED}");

  fprintf(outf, "\n\nthermo\t10000");
  fprintf(outf, "\nthermo_style\tcustom step time temp evdwl ecoul ebond");
  fprintf(outf, " eangle emol epair pe ke etotal fmax fnorm press");
  fprintf(outf, "\nthermo_modify\tflush yes");

  fprintf(outf, "\n\ntimestep\t%.6f", fib.parametersMD.timestep);
  fprintf(outf, "\nrun\t%.i upto", fib.parametersMD.runtime);

  fclose(outf);
}

void genQsub(collagenFibril fib)
{
  std::string file = filePaths.md_outputpath;
  file += "run.qsub";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "#!/bin/bash -l");
  fprintf(outf, "\n#$ -S /bin/bash");
  fprintf(outf, "\n#$ -l h_rt=%i:00:00", fib.parametersMD.walltime);
  fprintf(outf, "\n#$ -l mem=1G");
  // fprintf(outf, "\n#$ -N name"); // name is given automatically per script
  fprintf(outf, "\n#$ -pe mpi %i", fib.parametersMD.cores);
  fprintf(outf, "\n#$ -cwd");

  std::vector<md_var> md_vars = collectMDvars(fib);
  std::string dir = getDir(md_vars);

  fprintf(outf, "\n");
  fprintf(outf, "\ngerun %s", fib.parametersMD.lmp_mpi.c_str());
  fprintf(outf, " -in ./%s" , dir.c_str());
  fprintf(outf, "in.sim");

  fclose(outf);
}

void genBashScript(collagenFibril fib)
{
  std::string file = filePaths.md_outputpath;
  file += "sim.sh";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "#!/bin/bash -l");
  fprintf(outf, "\nLANG=en_US");

  fprintf(outf, "\n");
  fprintf(outf, "\nhelpFunction()");
  fprintf(outf, "\n{");
  fprintf(outf, "\n\techo \"\"");
  fprintf(outf, "\n\techo \"Usage: $0");
  fprintf(outf, " [-c]");
  fprintf(outf, " [-r]");
  fprintf(outf, "\"");
  fprintf(outf, "\n\techo -e \"\\t-c Creates LAMMPS output folders.");
  fprintf(outf, " Run this before option -r.\"");
  fprintf(outf, "\n\techo -e \"\\t-r Runs LAMMPS simulation on cluster.\"");
  fprintf(outf, "\n\texit 1");
  fprintf(outf, "\n}");

  fprintf(outf, "\n");
  fprintf(outf, "\ncreate=false");
  fprintf(outf, "\nrun=false");

  fprintf(outf, "\n");
  fprintf(outf, "\nwhile getopts \"cr\" opt");
  fprintf(outf, "\ndo");
  fprintf(outf, "\n\tcase \"$opt\" in");
  fprintf(outf, "\n\t\tc ) create=true ;;");
  fprintf(outf, "\n\t\tr ) run=true ;;");
  fprintf(outf, "\n\t\t? ) helpFunction ;;");
  fprintf(outf, "\n\tesac");
  fprintf(outf, "\ndone");


  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$create\" = false ]");
  fprintf(outf, " && [ \"$run\" = false ]");
  fprintf(outf, "\nthen");
  fprintf(outf, "\n\techo \"No option was chosen.\";");
  fprintf(outf, "\n\thelpFunction");
  fprintf(outf, "\nfi");

  std::vector<md_var> md_vars = collectMDvars(fib);
  std::string dir = getDir(md_vars);
  std::string cargs = getCargs(md_vars);
  std::string tabs = "\t";
  int i;

  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$create\" = true ]");
  fprintf(outf, "\nthen");
  for (i = 0; i < (int) md_vars.size(); i++) {
    tabs += "\t";
    printScriptVars(outf, md_vars[i], i + 1);
  }
  fprintf(outf, "\n%smkdir %s", tabs.c_str() , dir.c_str());
  fprintf(outf, "\n%scd %s", tabs.c_str(), filePaths.mainpath.c_str());
  fprintf(outf, "\n%s./main --config %s", tabs.c_str(), filePaths.configpath.c_str());
  fprintf(outf, " -md_sb");
  fprintf(outf, " -mdo %s -md_li %s", (filePaths.md_outputpath + dir).c_str(), cargs.c_str());
  fprintf(outf, "\n%scd %s", tabs.c_str(), filePaths.md_outputpath.c_str());
  for (int j = i; j >= 1; j--) {
    tabs.erase(0, 1);
    fprintf(outf, "\n%sdone", tabs.c_str());
  }
  fprintf(outf, "\nfi");

  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$run\" = true ]");
  fprintf(outf, "\nthen");
  for (i = 0; i < (int) md_vars.size(); i++) {
    tabs += "\t";
    printScriptVars(outf, md_vars[i], i + 1);
  }
  fprintf(outf, "\n%sqsub", tabs.c_str());
  std::replace(dir.begin(), dir.end(), '_', ',');
  if (!dir.empty()) {
    dir.pop_back();
    fprintf(outf, " -N %s", dir.c_str());
    fprintf(outf, " -v %s", dir.c_str());
  } else {
    fprintf(outf, " -N single_sim");
  }
  fprintf(outf, " run.qsub");
  for (int j = i; j >= 1; j--) {
    tabs.erase(0, 1);
    fprintf(outf, "\n%sdone", tabs.c_str());
  }
  fprintf(outf, "\nfi");

  fclose(outf);

  // make script executable
  chmod(file.c_str(), S_IRWXU);
}

void createLAMMPSfiles(collagenFibril fib)
{
  bool tmp = fib.parametersMD.outputTopology;
  tmp = tmp || fib.parametersMD.outputLAMMPSinput;
  tmp = tmp || fib.parametersMD.outputScript;

  if (flags.consoleOutput && tmp) {
    std::cout << "\n#\n#";
    std::cout << "\n# Writing LAMMPS files into folder ";
    std::cout << filePaths.md_outputpath;
  }

  if (fib.parametersMD.outputTopology) {
    if (fib.parametersMD.scriptbuild) {
      std::cout << "\n#  WARNING: Topology file is created multiple times!";
      std::cout << " To save disk space, it is recommended to avoid this by";
      std::cout << " unsetting MD_topology flag in the config file before";
      std::cout << " running sim.sh.";
    }
    if (flags.consoleOutput) {
      std::cout << "\n#  -> creating topology file for ";
      std::cout << pow(fib.parametersMD.numMolperDim, 3) << " molecules";
    }
    genTopologyZero(fib, fib.parametersMD.numMolperDim);
  }

  if (fib.parametersMD.outputLAMMPSinput) {
    if (flags.consoleOutput) {
      std::cout << "\n#  -> creating input file";
    }
    if (fib.parametersMD.scriptbuild) {
      genInSim(fib);
    } else {
      // TODO cycle through all variable parameters and create multiple in.sim
      // std::vector<md_var> md_vars = collectMDvars(fib);
      std::cout << " using *_start values for variable parameters";
      fib.parametersMD.kAngle = fib.parametersMD.kAngle_start;
      fib.parametersMD.dielectric = fib.parametersMD.dielectric_start;
      fib.parametersMD.LJepsilon = fib.parametersMD.LJepsilon_start;
      genInSim(fib);
    }
  }

  if (fib.parametersMD.outputScript) {
    if (fib.parametersMD.scriptbuild) {
      std::cout << "\n#  WARNING: Bash scripts are created multiple times!";
      std::cout << " To save disk space, it is recommended to avoid this by";
      std::cout << " unsetting MD_script flag in the config file before";
      std::cout << " running sim.sh.";
    }
    if (flags.consoleOutput) {
      std::cout << "\n#  -> creating bash scripts";
    }
    genBashScript(fib);
    genQsub(fib);
  }
}
