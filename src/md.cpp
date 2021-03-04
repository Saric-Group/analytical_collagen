#include "md.hpp"


/* Variables */
extern std::random_device device;
extern std::mt19937 generator;
extern filePaths_ filePaths;


/* Functions */
void genTopologyZero(collagenFibril fib, int L)
{
  // number of molecules per sidelength of cubic box
  // #define L 10

  // total number of molecules in box
  int numMol = L * L * L;
  // int numMol = 2;
  // std::cout << "\nnumMol=" << numMol;

  // defines the boundaries of cubic simulation box with length=boxlength
  // to be in line with Anne's volume fraction, we have to calculate:
  double result = 0.0;
  result = (36 - 1) * 0.255 * 1000 / (100 * 100 * 100);
  result /= (fib.mol.numAtoms - 1) * fib.mol.distanceAtoms * numMol;
  result = 1 / result;
  result = pow(result, 1. / 3.);
	double boxlength = result + fib.mol.length; //L * fib.mol.length;
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
  // std::cout << "\n Warning: boxlength=" << boxlength;
  // std::cout << "\n Warning: cube.length=" << grid.length;
  // std::cout << "\n Warning: cube.numPoints=" << grid.numPoints;
  fprintf(outf, "\n\nAtoms");
  // std::cout << "\ngridPoints=" << grid.gridPoints.size();
  std::uniform_real_distribution<double> dis_real(0.0, M_PI);
  double phi, theta;
  // std::cout << "\nmoleculeLength=" << fib.mol.length;
  fprintf(outf, "\n");
  for (int i = 0; i < numMol; i++) {
    phi = 0.;
    theta = 0.3 * M_PI;
    // phi = 2 * dis_real(generator);
    // theta = dis_real(generator);
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

void genInSim(collagenFibril fib)
{
  std::string file = filePaths.md_outputpath;
  file += "in.sim";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "log\tlog.sim");

  fprintf(outf, "\n\nunits\tlj");
  fprintf(outf,"\natom_style\tfull");

  fprintf(outf, "\n\nread_data\t");
  if (fib.parametersMD.scriptbuild) {
    fprintf(outf, "../");
  } else {
    std::string tmp = filePaths.md_outputpath;
    tmp.erase(0, 1);
    fprintf(outf, "/home/ucapkkl/Scratch/collagen");
    fprintf(outf, "%s", tmp.c_str());
  }
  fprintf(outf, "topology.0time");

  fprintf(outf, "\n\nvariable\tLANGEVIN_SEED\tequal\t1337");

  fprintf(outf, "\n\nbond_style\tharmonic");
  fprintf(outf, "\nbond_coeff\t1 500 %.6f", fib.mol.distanceAtoms);

  fprintf(outf, "\n\nangle_style\tharmonic");
  fprintf(outf, "\nangle_coeff\t1 50 180");

  fprintf(outf, "\n\npair_style\thybrid/overlay lj/cut %.3f coul/debye 1 %.3f", fib.parametersMD.lj_cutoff, fib.parametersMD.cd_cutoff);
  fprintf(outf, "\ndielectric\t%.3f", fib.parametersMD.dielectric);
  fprintf(outf, "\npair_coeff\t* * lj/cut %.3f 1 %.3f", fib.parametersMD.LJepsilon, fib.parametersMD.lj_cutoff);
  fprintf(outf, "\npair_coeff\t* * coul/debye %.3f", fib.parametersMD.cd_cutoff);
  fprintf(outf, "\npair_modify\tshift yes");
  fprintf(outf, "\nneigh_modify\texclude molecule/intra all");

  fprintf(outf, "\n\nvelocity\tall create 1.0 ${LANGEVIN_SEED}");
  fprintf(outf, "\nvelocity\tall zero angular");
  fprintf(outf, "\nvelocity\tall zero linear");

  fprintf(outf, "\n\ndump\tmydump all custom 10000 out_sim.xyz type x y z");
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
  fprintf(outf, "\n#$ -N");
  if (fib.parametersMD.scriptbuild) {
    fprintf(outf, " lje=%.3f_die=%.3f", fib.parametersMD.LJepsilon, fib.parametersMD.dielectric);
  } else {
    fprintf(outf, " Collagen_LAAMPS");
  }
  fprintf(outf, "\n#$ -pe mpi %i", fib.parametersMD.cores);
  fprintf(outf, "\n#$ -wd");
  std::string tmp = filePaths.md_outputpath;
  tmp.erase(0, 1);
  fprintf(outf, " /home/ucapkkl/Scratch/collagen");
  fprintf(outf, "%s", tmp.c_str());

  fprintf(outf, "\n\n");

  fprintf(outf, "\ngerun /home/ucapkkl/Scratch/lammps-29Oct20/src/lmp_mpi");
  fprintf(outf, " -in ./in.sim");
  fclose(outf);
}

void genBashScript()
{
  std::string file = filePaths.md_outputpath;
  file += "createFiles.sh";
  FILE *outf;
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "#!/bin/bash -l");
  fprintf(outf, "\nLANG=en_US");
  fprintf(outf, "\nfor lje in $(seq 0.01 0.1 0.21)");
  fprintf(outf, "\ndo");
  fprintf(outf, "\n\tfor die in $(seq 10.0 10.0 20.0)");
  fprintf(outf, "\n\tdo");
  fprintf(outf, "\n\t\tmkdir lje=${lje}_die=${die}");
  fprintf(outf, "\n\t\tcd ../../");
  fprintf(outf, "\n\t\t./main --config ./default.config");
  fprintf(outf, " -md_lje ${lje}");
  fprintf(outf, " -md_die ${die}");
  fprintf(outf, " -mdo %slje=${lje}_die=${die}/", filePaths.md_outputpath.c_str());
  fprintf(outf, " -md_li");
  fprintf(outf, " -md_sb");
  fprintf(outf, "\n\t\tcd %s", filePaths.md_outputpath.c_str());
  fprintf(outf, "\n\tdone");
  fprintf(outf, "\ndone");

  fclose(outf);

  // make script executable
  chmod(file.c_str(), S_IRWXU);


  file = filePaths.md_outputpath;
  file += "run.sh";
  outf = fopen(file.c_str(), "w");

  fprintf(outf, "#!/bin/bash -l");
  fprintf(outf, "\nLANG=en_US");
  fprintf(outf, "\nfor lje in $(seq 0.01 0.1 0.21)");
  fprintf(outf, "\ndo");
  fprintf(outf, "\n\tfor die in $(seq 10.0 10.0 20.0)");
  fprintf(outf, "\n\tdo");
  fprintf(outf, "\n\t\tcd ./lje=${lje}_die=${die}/");
  fprintf(outf, "\n\t\tqsub run.qsub");
  fprintf(outf, "\n\t\tcd ../");
  fprintf(outf, "\n\tdone");
  fprintf(outf, "\ndone");

  fclose(outf);

  // make script executable
  chmod(file.c_str(), S_IRWXU);
}

void createLAMMPSfiles(collagenFibril fib)
{
  if (fib.parametersMD.outputTopology) {
    if (fib.parametersMD.scriptbuild) {
      std::cout << "\n#\n#  WARNING: Topology file is created multiple times!";
      std::cout << " To save disk space, it is recommended to avoid this by";
      std::cout << " unsetting MD_topology flag in the config file before";
      std::cout << " running createFiles.sh.";
    }
    genTopologyZero(fib, fib.parametersMD.numMolperDim);
  }
  if (fib.parametersMD.outputLAMMPSinput) {
    genInSim(fib);
    genQsub(fib);
  }
  if (fib.parametersMD.outputScript) {
    if (fib.parametersMD.scriptbuild) {
      std::cout << "\n#\n#  WARNING: Bash scripts are created multiple times!";
      std::cout << " To save disk space, it is recommended to avoid this by";
      std::cout << " unsetting MD_script flag in the config file before";
      std::cout << " running createFiles.sh.";
    }
    genBashScript();
  }
}
