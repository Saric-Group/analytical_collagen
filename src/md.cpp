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
	double boxlength = 100.; //L * fib.mol.length;
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

  fprintf(outf, "\n\n%.3f %.3f xlo xhi", xlo, xhi);
  fprintf(outf, "\n%.3f %.3f ylo yhi", ylo, yhi);
  fprintf(outf, "\n%.3f %.3f zlo zhi", zlo, zhi);


  /* Masses */
  fprintf(outf, "\n\nMasses");
  fprintf(outf, "\n");
  for (int i = 0; i < fib.mol.numTypes; i++) {
    fprintf(outf, "\n\t%i 1.000", i + 1);
  }


  /* Atoms */
  cubeGrid grid(boxlength - fib.mol.length, L);
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
      if (fib.mol.atomTypes[j] == 21) {
        std::cout << "\n21 at " << i * fib.mol.numAtoms + j + 1;
      }
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

  fprintf(outf, "\n\nread_data\t/home/ucapkkl/Scratch/collagen");
  fprintf(outf, "%s", filePaths.md_outputpath.erase(0, 1).c_str());
  fprintf(outf, "topology.0time");

  fprintf(outf, "\n\nvariable\tLANGEVIN_SEED\tequal\t1337");

  fprintf(outf, "\n\nbond_style\tharmonic");
  fprintf(outf, "\nbond_coeff\t1 500 %.6f", fib.mol.distanceAtoms);

  fprintf(outf, "\n\nangle_style\tharmonic");
  fprintf(outf, "\nangle_coeff\t1 50 180");

  fprintf(outf, "\n\npair_style\thybrid/overlay lj/cut 5 coul/debye 1 5");
  fprintf(outf, "\ndielectric\t10");
  fprintf(outf, "\npair_coeff\t* * lj/cut 0.01 1 5");
  fprintf(outf, "\npair_coeff\t* * coul/debye 5");
  fprintf(outf, "\npair_modify\tshift yes");
  fprintf(outf, "\nneigh_modify\texclude molecule/intra all");

  fprintf(outf, "\n\nvelocity\tall create 1.0 ${LANGEVIN_SEED}");
  fprintf(outf, "\nvelocity\tall zero angular");
  fprintf(outf, "\nvelocity\tall zero linear");

  fprintf(outf, "\n\ndump\tmydump all custom 10000 out_sim.xyz type x y z");
  fprintf(outf, "\ndump_modify\tmydump sort id");

  // fprintf(outf, "\n\nfix\t0 all rigid/nve molecule");
  fprintf(outf, "\n\nfix\t0 all nve");
  fprintf(outf, "\nfix\t2 all langevin 1 1 1 ${LANGEVIN_SEED}");

  fprintf(outf, "\n\nthermo\t10000");
  fprintf(outf, "\nthermo_style\tcustom step time temp evdwl ecoul ebond");
  fprintf(outf, " eangle emol epair pe ke etotal fmax fnorm press");
  fprintf(outf, "\nthermo_modify\tflush yes");

  fprintf(outf, "\n\ntimestep\t0.002");
  fprintf(outf, "\nrun\t6000001 upto");

  fclose(outf);
}
