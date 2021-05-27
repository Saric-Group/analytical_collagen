#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "layermodel.hpp"
#include "random.hpp"
#include "funcs.hpp"
#include "md.hpp"
#include "cg.hpp"
#include "parse_config.hpp"



/* Settings */
filePaths_ filePaths;       /* IO File paths */
flags_ flags;               /* Various flags, see main.hpp */
dev_ dev;                   /* Variables used for dev stuff */



int main(int argc, char const *argv[])
{
  prioParse(argc, argv);
  programInfo();

  layerModel lm;
  collagenMolecule mol;

  if (parse(argc, argv, mol, lm) == 1) {
    return 1;
  }

  time_point start, end;
  getTime(start);

  printOptions(mol, lm);

  readAtomInfos(mol);
  lm.mol = mol;

  if (!flags.development) {
    // Layer-Model
    if (flags.minimize) {
      lm.minimizeEnergy();
      lm.writeConfig();
    }

    // LAMMPS
    if (mol.parametersMD.topology || mol.parametersMD.input) {
      // std::vector<double> coarseGrainedCharges;
      // coarseGrainedCharges = mol.smoothen(25);
      // coarseGrainedCharges = mol.binNormalize(coarseGrainedCharges, 100);
      // mol.readCharges(coarseGrainedCharges);
      // mol.printMoleculeInfo();
      createLAMMPSfiles(mol);
    }
  }

  // if (flags.md) {
  //   createLAMMPSfiles(mol);
  // }
  //
  // if (flags.layermodel) {
  //   lm.singleEmin();
  // }

  if (flags.development) {
    // mol.printMoleculeInfo();
    // lm.mol.numAtoms *= dev.factor;
    // lm.mol.numPos *= dev.factor;
    // lm.mol.numNeg *= dev.factor;
    // lm.minimizeEnergy();
    // lm.coutConfig();
    // std::string fdens = "./data/analysis/density";
    // lm.densityToFile(fdens);
    // lm.mol.numPos = dev.pos;
    // lm.mol.numNeg = dev.neg;
    // runRandomAnalysis(dev.samples, lm);
  }

  /************************************************************************/

  getTime(end);
  printCompTime(start, end);

  /************************************************************************/

  std::cout << "\n#\n# exiting...\n#\n";

  return 0;

}
