#include "parse_config.hpp"

using namespace config4cpp;


/* Variables */
extern filePaths_ filePaths;
extern flags_ flags;
extern dev_ dev;


/* Functions */
int parse(int argc, char const *argv[], collagenMolecule &mol, layerModel &lm)
{
  collagenMolecule dummyMol;
  layerModel dummyLm;

  std::cout << "\n# Parsing config file and command line arguments";
  for (int i = 0; i < argc; i++) {
    // scope deactivated for now
    // if (std::string(argv[i]) == "-scope") {
    //   if (i + 1 >= argc) {
    //     std::cerr << "Too few arguments after '-scope'" << "\n";
    //     return 1;
    //   }
    //   scope = argv[i + 1];
    // }
    if (std::string(argv[i]) == "-config") {
      if (i + 1 >= argc) {
        std::cerr << "\n#\n# Too few arguments after '-config'. Exit." << "\n";
        return 1;
      }
      filePaths.configpath = std::string(argv[i + 1]);
    }
  }

  setlocale(LC_ALL, "en_US.UTF-8");
  const char *scope = "collagen";
  // char scope[200];
  Configuration *cfg = Configuration::create();
  // SchemaValidator sv;

  for (int i = 0; i < argc; i++) {
    if (std::string(argv[i]) == "-set") {
      if (i + 2 >= argc) {
        std::cerr << "Too few arguments after '-set'" << "\n";
        return 1;
      }
      cfg->insertString(scope, argv[i + 1], argv[i + 2]);
    }
  }

  try {
    if (filePaths.configpath.empty()) {
      std::cout << "\n#  -> no config file selected, using defaults";
    } else {
      std::cout << "\n#  -> reading from config path " << filePaths.configpath;
      cfg->parse(filePaths.configpath.c_str());
      std::cout << "\n#  -> Options:";
    }
    // cfg->setFallbackConfiguration(Configuration::INPUT_STRING,
                                 // FallbackConfig::getString(), "default");
    // sv.parseSchema(FallbackConfig::getSchema());
    // sv.validate(cfg, "", "");
    // std::cout << "\n\n" << FallbackConfig::getString();

    // StringBuffer str;
    // cfg->dump(str, false);
    // std::cout << str.length();
    // the below prints all parsed variables
    // test = 1;
    // StringVector names;
    // cfg->listFullyScopedNames("", "", Configuration::CFG_SCOPE_AND_VARS, true, names);
    // for (int i = 0; i < names.length(); i++) {
    //   std::cout << "\n -> string " << i << ": " << names[i];
    // }
    // std::cout << "\n\n\nend";
    // std::cout.flush();
    // test = 2;
    /* General */
    // strcpy(scope, scope);
    // strcat(scope, ".general");
    // filePaths.mainpath = cfg->lookupString(scope, "general.mainpath");

    /* Flags */
    // strcpy(scope, scope);
    // strcat(scope, ".flags");
    flags.consoleOutput = cfg->lookupBoolean(scope, "flags.consoleOutput", false);
    flags.debugInfo = cfg->lookupBoolean(scope, "flags.debugInfo", false);
    flags.help = cfg->lookupBoolean(scope, "flags.help", false);
    flags.measureTime = cfg->lookupBoolean(scope, "flags.measureTime", false);

    /* io */
    // strcpy(scope, scope);
    // strcat(scope, ".io");
    filePaths.inputpath = cfg->lookupString(scope, "io.input", "./molecule/charges_distribution_36");
    filePaths.outputpath = cfg->lookupString(scope, "io.output", "./energy_min.dat");
    filePaths.md_outputpath = cfg->lookupString(scope, "io.lammps-output", "./md/");
    /* molecule */
    // strcpy(scope, scope);
    // strcat(scope, ".molecule");
    flags.printAtomInfo = cfg->lookupBoolean(scope, "molecule.printAtomInfo", false);
    flags.printMoleculeInfo = cfg->lookupBoolean(scope, "molecule.printMoleculeInfo", false);
    mol.diameterAtom = cfg->lookupFloat(scope, "molecule.diameterAtom", dummyMol.diameterAtom);
    mol.distanceAtoms = cfg->lookupFloat(scope, "molecule.distanceAtoms", dummyMol.distanceAtoms);

    /* main */
    flags.minimize = cfg->lookupBoolean(scope, "main.minimizeEnergy", false);
    lm.layers = cfg->lookupInt(scope, "main.layers", 2);
    lm.parameters.cd = cfg->lookupFloat(scope, "main.cd", dummyLm.parameters.cd);
    lm.parameters.cd_cutoff = cfg->lookupFloat(scope, "main.cdcutoff", dummyLm.parameters.cd_cutoff);
    lm.parameters.lj = cfg->lookupFloat(scope, "main.lj", dummyLm.parameters.lj);
    lm.parameters.lj_cutoff = cfg->lookupFloat(scope, "main.ljcutoff", dummyLm.parameters.lj_cutoff);
    lm.gap_stepsize = cfg->lookupFloat(scope, "main.radGapStepsize", dummyLm.gap_stepsize);
    lm.offset_stepsize = cfg->lookupFloat(scope, "main.offsetStepsize", dummyLm.offset_stepsize);

    // mol.parametersMD.script = cfg->lookupBoolean(scope, "md.script");
    mol.parametersMD.topology = cfg->lookupBoolean(scope, "main.lammps-topology", dummyMol.parametersMD.topology);
    mol.parametersMD.topFile = cfg->lookupString(scope, "main.topologyFile", dummyMol.parametersMD.topFile.c_str());

    mol.parametersMD.input = cfg->lookupBoolean(scope, "main.lammps-input", dummyMol.parametersMD.input);
    mol.parametersMD.inputFile = cfg->lookupString(scope, "main.inputFile", dummyMol.parametersMD.inputFile.c_str());
    mol.parametersMD.logFile = cfg->lookupString(scope, "main.logFile", dummyMol.parametersMD.logFile.c_str());
    mol.parametersMD.dumpFile = cfg->lookupString(scope, "main.dumpFile", dummyMol.parametersMD.dumpFile.c_str());

    mol.parametersMD.rigid = cfg->lookupBoolean(scope, "main.rigidMolecules", dummyMol.parametersMD.rigid);
    mol.parametersMD.random = cfg->lookupBoolean(scope, "main.randomOrientations", dummyMol.parametersMD.random);
    mol.parametersMD.phi = cfg->lookupFloat(scope, "main.phi", dummyMol.parametersMD.phi);
    mol.parametersMD.theta = cfg->lookupFloat(scope, "main.theta", dummyMol.parametersMD.theta);
    mol.parametersMD.numMolperDim = cfg->lookupInt(scope, "main.numMolperDim", dummyMol.parametersMD.numMolperDim);
    mol.parametersMD.timestep = cfg->lookupFloat(scope, "main.lammps-timestep", dummyMol.parametersMD.timestep);
    mol.parametersMD.runtime = cfg->lookupFloat(scope, "main.lammps-runtime", dummyMol.parametersMD.runtime);
    // strcat(scope, ".kAngle");
    mol.parametersMD.kAngle_start = cfg->lookupFloat(scope, "main.lammps-kAngle.start", dummyMol.parametersMD.kAngle_start);
    mol.parametersMD.kAngle_inc = cfg->lookupFloat(scope, "main.lammps-kAngle.increment", dummyMol.parametersMD.kAngle_inc);
    mol.parametersMD.kAngle_end = cfg->lookupFloat(scope, "main.lammps-kAngle.end", dummyMol.parametersMD.kAngle_end);
    // strcpy(scope, scope);
    // strcat(scope, ".MD");
    // strcat(scope, ".dielectric");
    mol.parametersMD.dielectric_start = cfg->lookupFloat(scope, "main.lammps-dielectric.start", dummyMol.parametersMD.dielectric_start);
    mol.parametersMD.dielectric_inc = cfg->lookupFloat(scope, "main.lammps-dielectric.increment", dummyMol.parametersMD.dielectric_inc);
    mol.parametersMD.dielectric_end = cfg->lookupFloat(scope, "main.lammps-dielectric.end", dummyMol.parametersMD.dielectric_end);
    mol.parametersMD.cd_cutoff = cfg->lookupFloat(scope, "main.lammps-dielectric.cutoff", dummyMol.parametersMD.cd_cutoff);
    // strcpy(scope, scope);
    // strcat(scope, ".MD");
    // strcat(scope, ".LJepsilon");
    mol.parametersMD.LJepsilon_start = cfg->lookupFloat(scope, "main.lammps-lj_epsilon.start", dummyMol.parametersMD.LJepsilon_start);
    mol.parametersMD.LJepsilon_inc = cfg->lookupFloat(scope, "main.lammps-lj_epsilon.increment", dummyMol.parametersMD.LJepsilon_inc);
    mol.parametersMD.LJepsilon_end = cfg->lookupFloat(scope, "main.lammps-lj_epsilon.end", dummyMol.parametersMD.LJepsilon_end);
    mol.parametersMD.lj_cutoff = cfg->lookupFloat(scope, "main.lammps-lj_epsilon.cutoff", dummyMol.parametersMD.lj_cutoff);

    /* layerModel */
    // // strcpy(scope, scope);
    // // strcat(scope, ".layerModel");
    // flags.layermodel = cfg->lookupBoolean(scope, "layerModel.active");
    // lm.parameters.chargehash = cfg->lookupBoolean(scope,
    //                                               "layerModel.chargehash");
    // lm.parameters.csv_output = cfg->lookupBoolean(scope, "layerModel.csv");
    // lm.parameters.xyz_outputs = cfg->lookupBoolean(scope, "layerModel.xyz");
    // lm.layers = cfg->lookupInt(scope, "layerModel.layers");
    // lm.gap_stepsize = cfg->lookupFloat(scope, "layerModel.radialGapStepsize");
    // lm.offset_stepsize = cfg->lookupFloat(scope, "layerModel.offsetStepsize");
    // // strcat(scope, ".LJparameters");
    // lm.parameters.lj_active = cfg->lookupBoolean(scope, "layerModel.lj_parameters.activate");
    // lm.parameters.lj_specific = cfg->lookupBoolean(scope, "layerModel.lj_parameters.specific");
    // lm.parameters.lj_stepsize = cfg->lookupFloat(scope, "layerModel.lj_parameters.stepsize");
    // lm.parameters.lj_steps = cfg->lookupInt(scope, "layerModel.lj_parameters.steps");
    // lm.parameters.lj_min = cfg->lookupFloat(scope, "layerModel.lj_parameters.min");
    // lm.parameters.lj_cutoff = cfg->lookupFloat(scope, "layerModel.lj_parameters.cutoff");
    // // strcpy(scope, scope);
    // // strcat(scope, ".layerModel");
    // // strcat(scope, ".CDparameters");
    // lm.parameters.cd_active = cfg->lookupBoolean(scope, "layerModel.cdParameters.activate");
    // lm.parameters.cd_specific = cfg->lookupBoolean(scope, "layerModel.cdParameters.specific");
    // lm.parameters.cd_stepsize = cfg->lookupFloat(scope, "layerModel.cdParameters.stepsize");
    // lm.parameters.cd_steps = cfg->lookupInt(scope, "layerModel.cdParameters.steps");
    // lm.parameters.cd_min = cfg->lookupFloat(scope, "layerModel.cdParameters.min");
    // lm.parameters.cd_cutoff = cfg->lookupFloat(scope, "layerModel.cdParameters.cutoff");

    /* MD */
    // strcpy(scope, scope);
    // strcat(scope, ".MD");
    // flags.md = cfg->lookupBoolean(scope, "md.active");
    // mol.parametersMD.lmp_mpi = cfg->lookupString(scope, "md.lmp_mpi");
    // mol.parametersMD.script = cfg->lookupBoolean(scope, "md.script");
    // mol.parametersMD.topology = cfg->lookupBoolean(scope, "md.topology");
    // mol.parametersMD.rigid = cfg->lookupBoolean(scope, "md.rigidMolecules");
    // mol.parametersMD.random = cfg->lookupBoolean(scope, "md.randomOrientations");
    // mol.parametersMD.phi = cfg->lookupFloat(scope, "md.phi");
    // mol.parametersMD.theta = cfg->lookupFloat(scope, "md.theta");
    // mol.parametersMD.input = cfg->lookupBoolean(scope, "md.input");
    // mol.parametersMD.walltime = cfg->lookupInt(scope, "md.walltime");
    // mol.parametersMD.cores = cfg->lookupInt(scope, "md.cores");
    // mol.parametersMD.numMolperDim = cfg->lookupInt(scope, "md.numMolperDim");
    // mol.parametersMD.timestep = cfg->lookupFloat(scope, "md.lammps_timestep");
    // mol.parametersMD.runtime = cfg->lookupFloat(scope, "md.lammps_runtime");
    // // strcat(scope, ".kAngle");
    // mol.parametersMD.kAngle_start = cfg->lookupFloat(scope, "md.kAngle.start");
    // mol.parametersMD.kAngle_inc = cfg->lookupFloat(scope, "md.kAngle.increment");
    // mol.parametersMD.kAngle_end = cfg->lookupFloat(scope, "md.kAngle.end");
    // // strcpy(scope, scope);
    // // strcat(scope, ".MD");
    // // strcat(scope, ".dielectric");
    // mol.parametersMD.dielectric_start = cfg->lookupFloat(scope, "md.dielectric.start");
    // mol.parametersMD.dielectric_inc = cfg->lookupFloat(scope, "md.dielectric.increment");
    // mol.parametersMD.dielectric_end = cfg->lookupFloat(scope, "md.dielectric.end");
    // mol.parametersMD.cd_cutoff = cfg->lookupFloat(scope, "md.dielectric.cutoff");
    // // strcpy(scope, scope);
    // // strcat(scope, ".MD");
    // // strcat(scope, ".LJepsilon");
    // mol.parametersMD.LJepsilon_start = cfg->lookupFloat(scope, "md.lj_epsilon.start");
    // mol.parametersMD.LJepsilon_inc = cfg->lookupFloat(scope, "md.lj_epsilon.increment");
    // mol.parametersMD.LJepsilon_end = cfg->lookupFloat(scope, "md.lj_epsilon.end");
    // mol.parametersMD.lj_cutoff = cfg->lookupFloat(scope, "md.lj_epsilon.cutoff");

    /* development */
    // strcpy(scope, scope);
    // strcat(scope, ".development");
    flags.development = cfg->lookupBoolean(scope, "development.active", false);
    dev.random = cfg->lookupBoolean(scope, "development.random", false);
    dev.seed = cfg->lookupInt(scope, "development.seed", 1337);
    dev.factor = cfg->lookupFloat(scope, "development.factor", 1.0);
    dev.samples = cfg->lookupInt(scope, "development.samples", 1);
    dev.pos = cfg->lookupInt(scope, "development.pos", 86);
    dev.neg = cfg->lookupInt(scope, "development.neg", 82);

  } catch (const ConfigurationException &ex) {
    std::cerr << "\n# error: " << ex.c_str() << "\n";
    cfg->destroy();
  }

  cfg->destroy();

  return 0;
}
