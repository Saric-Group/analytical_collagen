#include "md.hpp"


/* Variables */
extern std::random_device device;
extern std::mt19937 generator;
extern filePaths_ filePaths;
extern flags_ flags;


/* Functions */
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
  // Increased the box, lowered the volume fraction to keep molecules in box
  // might need to think about boxsize again
	double boxlength = result + fib.mol.length;
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
  fprintf(outf, "\n\nAtoms");
  std::uniform_real_distribution<double> dis_real(0.0, M_PI);
  double phi, theta;
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

  fprintf(outf, "log\t");
  fprintf(outf, "./lje=%.3f_die=%.1f", fib.parametersMD.LJepsilon, fib.parametersMD.dielectric);
  fprintf(outf, "_ka=%.1f", fib.parametersMD.kAngle);
  fprintf(outf, "_cdc=%.1f_ljc=%.1f", fib.parametersMD.cd_cutoff, fib.parametersMD.lj_cutoff);
  fprintf(outf, "_ts=%.6f_rt=%.i", fib.parametersMD.timestep, fib.parametersMD.runtime);
  fprintf(outf, "/log.sim");

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

  fprintf(outf, "\n\ndump\tmydump all custom 10000");
  fprintf(outf, " ./lje=%.3f_die=%.1f", fib.parametersMD.LJepsilon, fib.parametersMD.dielectric);
  fprintf(outf, "_ka=%.1f", fib.parametersMD.kAngle);
  fprintf(outf, "_cdc=%.1f_ljc=%.1f", fib.parametersMD.cd_cutoff, fib.parametersMD.lj_cutoff);
  fprintf(outf, "_ts=%.6f_rt=%.i", fib.parametersMD.timestep, fib.parametersMD.runtime);
  fprintf(outf, "/out_sim.xyz type q x y z");
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

  fprintf(outf, "\n");
  // fprintf(outf, "\ngerun /home/ucapkkl/Scratch/lammps-29Oct20/src/lmp_mpi");
  fprintf(outf, "\ngerun %s", fib.parametersMD.lmp_mpi.c_str());
  fprintf(outf, " -in ./lje=${lje}_die=${die}");
  fprintf(outf, "_ka=${ka}");
  fprintf(outf, "_cdc=${cdc}_ljc=${ljc}");
  fprintf(outf, "_ts=${ts}_rt=${rt}/");
  fprintf(outf, "in.sim");

  fclose(outf);
}

void genBashScript(collagenFibril fib)
{
  // might want to use subfunctions for the creation of the script files
  // but for now this suffices
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
  // fprintf(outf, "\n\techo \"Usage: $0");
  // fprintf(outf, " -a paramA -b paramB -c paramC\"");
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
  // fprintf(outf, "\nwhile getopts \"a:b:c:\" opt");
  fprintf(outf, "\ndo");
  fprintf(outf, "\n\tcase \"$opt\" in");
  fprintf(outf, "\n\t\tc ) create=true ;;");
  fprintf(outf, "\n\t\tr ) run=true ;;");
  // fprintf(outf, "\n\t\ta ) paramA=\"$OPTARG\" ;;");
  // fprintf(outf, "\n\t\tb ) paramB=\"$OPTARG\" ;;");
  // fprintf(outf, "\n\t\tc ) paramC=\"$OPTARG\" ;;");
  fprintf(outf, "\n\t\t? ) helpFunction ;;");
  fprintf(outf, "\n\tesac");
  fprintf(outf, "\ndone");

  // fprintf(outf, "\n");
  // fprintf(outf, "\nif [ -z \"$paramA\" ] ");
  // fprintf(outf, "|| [ -z \"$paramB\" ]");
  // fprintf(outf, "|| [ -z \"$paramC\" ]");
  // fprintf(outf, "\nthen");
  // fprintf(outf, "\n\techo \"Some or all parameters are empty\";");
  // fprintf(outf, "\n\thelpFunction");
  // fprintf(outf, "\nfi");

  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$create\" = false ]");
  fprintf(outf, " && [ \"$run\" = false ]");
  fprintf(outf, "\nthen");
  fprintf(outf, "\n\techo \"No option was chosen.\";");
  fprintf(outf, "\n\thelpFunction");
  fprintf(outf, "\nfi");

  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$create\" = true ]");
  fprintf(outf, "\nthen");
  fprintf(outf, "\n\tfor lje in $(seq");
  fprintf(outf, " %.3f", fib.parametersMD.LJepsilon_start);
  fprintf(outf, " %.3f", fib.parametersMD.LJepsilon_inc);
  fprintf(outf, " %.3f)", fib.parametersMD.LJepsilon_end);
  fprintf(outf, "\n\tdo");
  fprintf(outf, "\n\t\tfor die in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.dielectric_start);
  fprintf(outf, " %.1f", fib.parametersMD.dielectric_inc);
  fprintf(outf, " %.1f)", fib.parametersMD.dielectric_end);
  fprintf(outf, "\n\t\tdo");
  fprintf(outf, "\n\t\t\tfor ka in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.kAngle_start);
  fprintf(outf, " %.1f", fib.parametersMD.kAngle_inc);
  fprintf(outf, " %.1f)", fib.parametersMD.kAngle_end);
  fprintf(outf, "\n\t\t\tdo");
  fprintf(outf, "\n\t\t\t\tfor cdc in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.cd_cutoff);
  fprintf(outf, " %.1f", fib.parametersMD.cd_cutoff);
  fprintf(outf, " %.1f)", fib.parametersMD.cd_cutoff);
  fprintf(outf, "\n\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\tfor ljc in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.lj_cutoff);
  fprintf(outf, " %.1f", fib.parametersMD.lj_cutoff);
  fprintf(outf, " %.1f)", fib.parametersMD.lj_cutoff);
  fprintf(outf, "\n\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\tfor ts in $(seq");
  fprintf(outf, " %.6f", fib.parametersMD.timestep);
  fprintf(outf, " %.6f", fib.parametersMD.timestep);
  fprintf(outf, " %.6f)", fib.parametersMD.timestep);
  fprintf(outf, "\n\t\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\t\tfor rt in $(seq");
  fprintf(outf, " %.i", fib.parametersMD.runtime);
  fprintf(outf, " %.i", fib.parametersMD.runtime);
  fprintf(outf, " %.i)", fib.parametersMD.runtime);
  fprintf(outf, "\n\t\t\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\t\t\tmkdir lje=${lje}_die=${die}");
  fprintf(outf, "_ka=${ka}");
  fprintf(outf, "_cdc=${cdc}_ljc=${ljc}");
  fprintf(outf, "_ts=${ts}_rt=${rt}");
  fprintf(outf, "\n\t\t\t\t\t\t\t\tcd ../../");
  fprintf(outf, "\n\t\t\t\t\t\t\t\t./main --config ./default.config");
  fprintf(outf, " -mdo %s", filePaths.md_outputpath.c_str());
  fprintf(outf, "lje=${lje}_die=${die}");
  fprintf(outf, "_ka=${ka}");
  fprintf(outf, "_cdc=${cdc}_ljc=${ljc}");
  fprintf(outf, "_ts=${ts}_rt=${rt}/");
  fprintf(outf, " -md_li");
  fprintf(outf, " -md_sb");
  fprintf(outf, " -md_lje ${lje}");
  fprintf(outf, " -md_die ${die}");
  fprintf(outf, " -md_ka ${ka}");
  fprintf(outf, " -md_cdc ${cdc}");
  fprintf(outf, " -md_ljc ${ljc}");
  fprintf(outf, " -md_ts ${ts}");
  fprintf(outf, " -md_rt ${rt}");
  fprintf(outf, "\n\t\t\t\t\t\t\t\tcd %s", filePaths.md_outputpath.c_str());
  fprintf(outf, "\n\t\t\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\tdone");
  fprintf(outf, "\n\t\t\tdone");
  fprintf(outf, "\n\t\tdone");
  fprintf(outf, "\n\tdone");
  fprintf(outf, "\nfi");

  fprintf(outf, "\n");
  fprintf(outf, "\nif [ \"$run\" = true ]");
  fprintf(outf, "\nthen");
  fprintf(outf, "\n\tfor lje in $(seq");
  fprintf(outf, " %.3f", fib.parametersMD.LJepsilon_start);
  fprintf(outf, " %.3f", fib.parametersMD.LJepsilon_inc);
  fprintf(outf, " %.3f)", fib.parametersMD.LJepsilon_end);
  fprintf(outf, "\n\tdo");
  fprintf(outf, "\n\t\tfor die in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.dielectric_start);
  fprintf(outf, " %.1f", fib.parametersMD.dielectric_inc);
  fprintf(outf, " %.1f)", fib.parametersMD.dielectric_end);
  fprintf(outf, "\n\t\tdo");
  fprintf(outf, "\n\t\t\tfor ka in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.kAngle_start);
  fprintf(outf, " %.1f", fib.parametersMD.kAngle_inc);
  fprintf(outf, " %.1f)", fib.parametersMD.kAngle_end);
  fprintf(outf, "\n\t\t\tdo");
  fprintf(outf, "\n\t\t\t\tfor cdc in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.cd_cutoff);
  fprintf(outf, " %.1f", fib.parametersMD.cd_cutoff);
  fprintf(outf, " %.1f)", fib.parametersMD.cd_cutoff);
  fprintf(outf, "\n\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\tfor ljc in $(seq");
  fprintf(outf, " %.1f", fib.parametersMD.lj_cutoff);
  fprintf(outf, " %.1f", fib.parametersMD.lj_cutoff);
  fprintf(outf, " %.1f)", fib.parametersMD.lj_cutoff);
  fprintf(outf, "\n\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\tfor ts in $(seq");
  fprintf(outf, " %.6f", fib.parametersMD.timestep);
  fprintf(outf, " %.6f", fib.parametersMD.timestep);
  fprintf(outf, " %.6f)", fib.parametersMD.timestep);
  fprintf(outf, "\n\t\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\t\tfor rt in $(seq");
  fprintf(outf, " %.i", fib.parametersMD.runtime);
  fprintf(outf, " %.i", fib.parametersMD.runtime);
  fprintf(outf, " %.i)", fib.parametersMD.runtime);
  fprintf(outf, "\n\t\t\t\t\t\t\tdo");
  fprintf(outf, "\n\t\t\t\t\t\t\t\tqsub");
  fprintf(outf, " -N l=${lje}_c=${die}");
  fprintf(outf, "_ka=${ka}");
  fprintf(outf, "_lc=${ljc}_cc=${cdc}");
  fprintf(outf, "_ts=${ts}_rt=${rt}");
  fprintf(outf, " -v lje=${lje},die=${die}");
  fprintf(outf, ",ka=${ka}");
  fprintf(outf, ",cdc=${cdc},ljc=${ljc}");
  fprintf(outf, ",ts=${ts},rt=${rt}");
  fprintf(outf, " run.qsub");
  fprintf(outf, "\n\t\t\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\t\tdone");
  fprintf(outf, "\n\t\t\t\tdone");
  fprintf(outf, "\n\t\t\tdone");
  fprintf(outf, "\n\t\tdone");
  fprintf(outf, "\n\tdone");
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
    genInSim(fib);
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
