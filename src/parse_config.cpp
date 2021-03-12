#include "parse_config.hpp"

using namespace config4cpp;


/* Variables */
extern filePaths_ filePaths;
extern flags_ flags;


/* Functions */
void parse_config(collagenMolecule mol)
{
  setlocale(LC_ALL, "en_US.UTF-8");
  std::string tmp = "";

  const char *scope = "collagen";
  Configuration *cfg = Configuration::create();

  try {
    cfg->parse("./config/test.config");

    /* General */
    filePaths.mainpath = cfg->lookupString(scope, "general.mainpath");

    /* Flags */
    flags.consoleOutput = cfg->lookupBoolean(scope, "flags.consoleOutput");
    flags.help = cfg->lookupBoolean(scope, "flags.help");
    flags.measureTime = cfg->lookupBoolean(scope, "flags.measureTime");

    /* io */
    filePaths.inputpath = cfg->lookupString(scope, "io.input");
    filePaths.outputpath = cfg->lookupString(scope, "io.output");
    filePaths.md_outputpath = cfg->lookupString(scope, "io.mdoutput");

    /* molecule */
    flags.printAtomInfo = cfg->lookupBoolean(scope, "molecule.printAtomInfo");
    tmp = "molecule.printMoleculeInfo";
    flags.printMoleculeInfo = cfg->lookupBoolean(scope, tmp.c_str());
    mol.diameterAtom = cfg->lookupFloat(scope, "molecule.diameterAtom");
    mol.distanceAtoms = cfg->lookupFloat(scope, "molecule.distanceAtoms");

    /* layerModel */

  } catch (const ConfigurationException &ex) {
    std::cerr << ex.c_str() << "\n";
    cfg->destroy();
  }

  cfg->destroy();
}
