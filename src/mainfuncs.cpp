#include "mainfuncs.hpp"


/* Variables */
overrides_ overrides;
extern filePaths_ filePaths;
extern flags_ flags;

/* Functions */
bool fexists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
      fclose(file);
      return true;
  } else {
      return false;
  }
}

int process_arg(std::string strarg, int* errcodes, bool dashes, layerModel &lm, collagenMolecule &mol)
{
  //io overrides
    int errstate = 0;

    std::string d = dashes ? "-" : "";
    std::string dd = dashes ? "--" : "";
    if(overrides.mainPath)
    {
      if(verify_path(strarg,errcodes) == 0) {
        filePaths.mainpath = strarg;
      }
      overrides.mainPath = false;
    }
    if(overrides.input)
    {
      if(verify_path(strarg,errcodes) == 0) {
        filePaths.inputpath = strarg;
      }
      overrides.input = false;
    }

    if(overrides.output)
    {
      filePaths.outputpath = strarg;
      overrides.output = false;
    }

    if(overrides.md_output)
    {
      filePaths.md_outputpath = strarg;
      overrides.md_output = false;
    }

    //spatial overrides

    if(overrides.layers)
    {
      lm.layers = safe_read_integer(lm.layers,strarg,errcodes);
      overrides.layers = false;
    }
    if(overrides.diameter)
    {
      lm.mol.diameterAtom = safe_read_double(lm.mol.diameterAtom,strarg,errcodes);
      overrides.diameter = false;
    }
    if(overrides.distance)
    {
      lm.mol.distanceAtoms = safe_read_double(lm.mol.distanceAtoms,strarg,errcodes);
      overrides.distance = false;
    }

    //potential overrides

    if(overrides.ljsteps)
    {
      lm.parameters.lj_steps = safe_read_integer(lm.parameters.lj_steps,strarg,errcodes);
      overrides.ljsteps = false;
    }
    if(overrides.ljstepsize)
    {
      lm.parameters.lj_stepsize = safe_read_double(lm.parameters.lj_stepsize,strarg,errcodes);
      overrides.ljstepsize = false;
    }
    if(overrides.cdsteps)
    {
      lm.parameters.cd_steps = safe_read_integer(lm.parameters.cd_steps,strarg,errcodes);
      overrides.cdsteps = false;
    }
    if(overrides.cdstepsize)
    {
      lm.parameters.cd_stepsize = safe_read_double(lm.parameters.cd_stepsize,strarg,errcodes);
      overrides.cdstepsize = false;
    }
    if(overrides.ljmin)
    {
      lm.parameters.lj_min = safe_read_double(lm.parameters.lj_min,strarg,errcodes);
      overrides.ljmin = false;
    }
    if(overrides.cdmin)
    {
      lm.parameters.cd_min = safe_read_double(lm.parameters.cd_min,strarg,errcodes);
      overrides.cdmin = false;
    }
    if(overrides.ljcut)
    {
      lm.parameters.lj_cutoff = safe_read_double(lm.parameters.lj_cutoff,strarg,errcodes);
      lm.parameters.calcCutoff();
      overrides.ljcut = false;
    }
    if(overrides.cdcut)
    {
      lm.parameters.cd_cutoff = safe_read_double(lm.parameters.cd_cutoff,strarg,errcodes);
      lm.parameters.calcCutoff();
      overrides.cdcut = false;
    }

    // MD related
    if(overrides.MD_mpi_dir)
    {
      mol.parametersMD.lmp_mpi = strarg;
      overrides.MD_mpi_dir = false;
    }
    if(overrides.MD_walltime)
    {
      mol.parametersMD.walltime = safe_read_integer(mol.parametersMD.walltime,strarg,errcodes);
      overrides.MD_walltime = false;
    }
    if(overrides.MD_cores)
    {
      mol.parametersMD.cores = safe_read_integer(mol.parametersMD.cores,strarg,errcodes);
      overrides.MD_cores = false;
    }
    if(overrides.MD_numMolperDim)
    {
      mol.parametersMD.numMolperDim = safe_read_integer(mol.parametersMD.numMolperDim,strarg,errcodes);
      overrides.MD_numMolperDim = false;
    }
    if(overrides.MD_phi)
    {
      mol.parametersMD.phi = safe_read_double(mol.parametersMD.phi,strarg,errcodes);
      overrides.MD_phi = false;
    }
    if(overrides.MD_theta)
    {
      mol.parametersMD.theta = safe_read_double(mol.parametersMD.theta,strarg,errcodes);
      overrides.MD_theta = false;
    }

    if(overrides.MD_kAngle)
    {
      mol.parametersMD.kAngle = safe_read_double(mol.parametersMD.kAngle,strarg,errcodes);
      overrides.MD_kAngle = false;
    }
    if(overrides.MD_kAngle_start)
    {
      mol.parametersMD.kAngle_start = safe_read_double(mol.parametersMD.kAngle_start,strarg,errcodes);
      overrides.MD_kAngle_start = false;
    }
    if(overrides.MD_kAngle_inc)
    {
      mol.parametersMD.kAngle_inc = safe_read_double(mol.parametersMD.kAngle_inc,strarg,errcodes);
      overrides.MD_kAngle_inc = false;
    }
    if(overrides.MD_kAngle_end)
    {
      mol.parametersMD.kAngle_end = safe_read_double(mol.parametersMD.kAngle_end,strarg,errcodes);
      overrides.MD_kAngle_end = false;
    }

    if(overrides.MD_dielectric)
    {
      mol.parametersMD.dielectric = safe_read_double(mol.parametersMD.dielectric,strarg,errcodes);
      overrides.MD_dielectric = false;
    }
    if(overrides.MD_dielectric_start)
    {
      mol.parametersMD.dielectric_start = safe_read_double(mol.parametersMD.dielectric_start,strarg,errcodes);
      overrides.MD_dielectric_start = false;
    }
    if(overrides.MD_dielectric_inc)
    {
      mol.parametersMD.dielectric_inc = safe_read_double(mol.parametersMD.dielectric_inc,strarg,errcodes);
      overrides.MD_dielectric_inc = false;
    }
    if(overrides.MD_dielectric_end)
    {
      mol.parametersMD.dielectric_end = safe_read_double(mol.parametersMD.dielectric_end,strarg,errcodes);
      overrides.MD_dielectric_end = false;
    }
    if(overrides.MD_LJepsilon)
    {
      mol.parametersMD.LJepsilon = safe_read_double(mol.parametersMD.LJepsilon,strarg,errcodes);
      overrides.MD_LJepsilon = false;
    }
    if(overrides.MD_LJepsilon_start)
    {
      mol.parametersMD.LJepsilon_start = safe_read_double(mol.parametersMD.LJepsilon_start,strarg,errcodes);
      overrides.MD_LJepsilon_start = false;
    }
    if(overrides.MD_LJepsilon_inc)
    {
      mol.parametersMD.LJepsilon_inc = safe_read_double(mol.parametersMD.LJepsilon_inc,strarg,errcodes);
      overrides.MD_LJepsilon_inc = false;
    }
    if(overrides.MD_LJepsilon_end)
    {
      mol.parametersMD.LJepsilon_end = safe_read_double(mol.parametersMD.LJepsilon_end,strarg,errcodes);
      overrides.MD_LJepsilon_end = false;
    }
    if(overrides.MD_LJcutoff)
    {
      mol.parametersMD.lj_cutoff = safe_read_double(mol.parametersMD.lj_cutoff,strarg,errcodes);
      overrides.MD_LJcutoff = false;
    }
    if(overrides.MD_CDcutoff)
    {
      mol.parametersMD.cd_cutoff = safe_read_double(mol.parametersMD.cd_cutoff,strarg,errcodes);
      overrides.MD_CDcutoff = false;
    }
    if(overrides.MD_timestep)
    {
      mol.parametersMD.timestep = safe_read_double(mol.parametersMD.timestep,strarg,errcodes);
      overrides.MD_timestep = false;
    }
    if(overrides.MD_runtime)
    {
      mol.parametersMD.runtime = safe_read_integer(mol.parametersMD.runtime,strarg,errcodes);
      overrides.MD_runtime = false;
    }

    //error handling

    if (parse_errs(*errcodes,strarg,true) != 0)
    {
      errstate = 1;
    }



    //bool flags and help

    if(flag(strarg,d+"co") || flag(strarg,dd+"consoleOutput"))
    {
      flags.consoleOutput = true;
    }

    if(flag(strarg,d+"h") || flag(strarg,dd+"help"))
    {
      print_help();
      return -1;
    }

    if(flag(strarg,d+"c") || flag(strarg,dd+"chargehash"))
    {
      lm.parameters.chargehash = true;
    }

    if(flag(strarg,d+"x") || flag(strarg,dd+"xyz"))
    {
      lm.parameters.xyz_outputs = true;
    }

    if(flag(strarg,d+"csv") || flag(strarg,dd+"csv"))
    {
      lm.parameters.csv_output = true;
    }

    if(flag(strarg,d+"pai") || flag(strarg,dd+"printAtomInfo"))
    {
      flags.printAtomInfo = true;
    }

    if(flag(strarg,d+"pmi") || flag(strarg,dd+"printMoleculeInfo"))
    {
      flags.printMoleculeInfo = true;
    }

    if(flag(strarg,d+"oo") || flag(strarg,dd+"originalOutput"))
    {
      flags.layermodel = true;
    }

    if(flag(strarg,d+"md") || flag(strarg,dd+"MD"))
    {
      flags.md = true;
    }

    if(flag(strarg,d+"time") || flag(strarg,dd+"measureTime"))
    {
      flags.measureTime = true;
    }

    if(flag(strarg,d+"dev") || flag(strarg,dd+"development"))
    {
      flags.development = true;
    }

    // general overrides
    if(flag(strarg,d+"mp") || flag(strarg,dd+"mainPath"))
    {
      overrides.mainPath = true;
    }

    //io overrides

    if(flag(strarg,d+"i") || flag(strarg,dd+"input"))
    {
      overrides.input = true;
    }

    if(flag(strarg,d+"o") || flag(strarg,dd+"output"))
    {
      overrides.output = true;
    }
    if(flag(strarg,d+"mdo") || flag(strarg,dd+"mdoutput"))
    {
      overrides.md_output = true;
    }

    //spatial overrides

    if(flag(strarg,d+"l") || flag(strarg,dd+"layers"))
    {
      overrides.layers = true;
    }
    if(flag(strarg,d+"dia") || flag(strarg,dd+"diameter"))
    {
      overrides.diameter = true;
    }
    if(flag(strarg,d+"d") || flag(strarg,dd+"dist"))
    {
      overrides.distance = true;
    }

    //potential overrides

    if(flag(strarg,d+"slj") || flag(strarg,dd+"ljsteps"))
    {
      overrides.ljsteps = true;
    }
    if(flag(strarg,d+"sslj") || flag(strarg,dd+"ljstepssize"))
    {
      overrides.ljstepsize = true;
    }
    if(flag(strarg,d+"scd") || flag(strarg,dd+"cdsteps"))
    {
      overrides.cdsteps = true;
    }
    if(flag(strarg,d+"sscd") || flag(strarg,dd+"cdstepssize"))
    {
      overrides.cdstepsize = true;
    }
    if(flag(strarg,d+"mcd") || flag(strarg,dd+"cdmin"))
    {
      overrides.cdmin = true;
    }
    if(flag(strarg,d+"mlj") || flag(strarg,dd+"ljmin"))
    {
      overrides.ljmin = true;
    }


    if(flag(strarg,d+"ccd") || flag(strarg,dd+"cdcut"))
    {
      overrides.cdcut = true;
    }
    if(flag(strarg,d+"clj") || flag(strarg,dd+"ljcut"))
    {
      overrides.ljcut = true;
    }

    // MD related
    if(flag(strarg,d+"md_s") || flag(strarg,dd+"MD_script"))
    {
      mol.parametersMD.script = true;
    }
    if(flag(strarg,d+"md_t") || flag(strarg,dd+"MD_topology"))
    {
      mol.parametersMD.topology = true;
    }
    if(flag(strarg,d+"md_li") || flag(strarg,dd+"MD_LAMMPS_input"))
    {
      mol.parametersMD.input = true;
    }
    if(flag(strarg,d+"md_sb"))
    {
      mol.parametersMD.scriptbuild = true;
    }
    if(flag(strarg,d+"md_rigid") || flag(strarg,dd+"MD_rigid"))
    {
      mol.parametersMD.rigid = true;
    }
    if(flag(strarg,d+"md_mpid") || flag(strarg,dd+"MD_lmp_mpi"))
    {
      overrides.MD_mpi_dir = true;
    }
    if(flag(strarg,d+"md_wt") || flag(strarg,dd+"MD_walltime"))
    {
      overrides.MD_walltime = true;
    }
    if(flag(strarg,d+"md_c") || flag(strarg,dd+"MD_cores"))
    {
      overrides.MD_cores = true;
    }
    if(flag(strarg,d+"md_N") || flag(strarg,dd+"MD_N"))
    {
      overrides.MD_numMolperDim = true;
    }
    if(flag(strarg,d+"md_ran") || flag(strarg,dd+"MD_random"))
    {
      mol.parametersMD.random = true;
    }
    if(flag(strarg,d+"md_phi") || flag(strarg,dd+"MD_phi"))
    {
      overrides.MD_phi = true;
    }
    if(flag(strarg,d+"md_the") || flag(strarg,dd+"MD_theta"))
    {
      overrides.MD_theta = true;
    }
    if(flag(strarg,d+"md_ka") || flag(strarg,dd+"MD_kAngle"))
    {
      overrides.MD_kAngle = true;
    }
    if(flag(strarg,d+"md_kas") || flag(strarg,dd+"MD_kAngle_start"))
    {
      overrides.MD_kAngle_start = true;
    }
    if(flag(strarg,d+"md_kai") || flag(strarg,dd+"MD_kAngle_inc"))
    {
      overrides.MD_kAngle_inc = true;
    }
    if(flag(strarg,d+"md_kae") || flag(strarg,dd+"MD_kAngle_end"))
    {
      overrides.MD_kAngle_end = true;
    }
    if(flag(strarg,d+"md_die") || flag(strarg,dd+"MD_dielectric"))
    {
      overrides.MD_dielectric = true;
    }
    if(flag(strarg,d+"md_ds") || flag(strarg,dd+"MD_dielectric_start"))
    {
      overrides.MD_dielectric_start = true;
    }
    if(flag(strarg,d+"md_di") || flag(strarg,dd+"MD_dielectric_inc"))
    {
      overrides.MD_dielectric_inc = true;
    }
    if(flag(strarg,d+"md_de") || flag(strarg,dd+"MD_dielectric_end"))
    {
      overrides.MD_dielectric_end = true;
    }
    if(flag(strarg,d+"md_lje") || flag(strarg,dd+"MD_LJepsilon"))
    {
      overrides.MD_LJepsilon = true;
    }
    if(flag(strarg,d+"md_es") || flag(strarg,dd+"MD_LJepsilon_start"))
    {
      overrides.MD_LJepsilon_start = true;
    }
    if(flag(strarg,d+"md_ei") || flag(strarg,dd+"MD_LJepsilon_inc"))
    {
      overrides.MD_LJepsilon_inc = true;
    }
    if(flag(strarg,d+"md_ee") || flag(strarg,dd+"MD_LJepsilon_end"))
    {
      overrides.MD_LJepsilon_end = true;
    }
    if(flag(strarg,d+"md_ljc") || flag(strarg,dd+"MD_LJcutoff"))
    {
      overrides.MD_LJcutoff = true;
    }
    if(flag(strarg,d+"md_cdc") || flag(strarg,dd+"MD_CDcutoff"))
    {
      overrides.MD_CDcutoff = true;
    }
    if(flag(strarg,d+"md_ts") || flag(strarg,dd+"MD_timestep"))
    {
      overrides.MD_timestep = true;
    }
    if(flag(strarg,d+"md_rt") || flag(strarg,dd+"MD_runtime"))
    {
      overrides.MD_runtime = true;
    }

    return errstate;
}

int read_config_file(std::string path, layerModel &lm, collagenMolecule &mol)
{
  int errcodes = 0;
  std::ifstream infile(path);
  std::string line;
  while (std::getline(infile, line))
  {
      errcodes = 0;
      std::string s = remove_spaces(line);
      if (s.length() > 0 && s.at(0) == '#') {
        continue;
      }
      std::string delimiter = "=";

      size_t pos = 0;
      std::string token;
      while ((pos = s.find(delimiter)) != std::string::npos) {
          token = s.substr(0, pos);
          process_arg(token, &errcodes, false, lm, mol);
          s.erase(0, pos + delimiter.length());
      }
      process_arg(s, &errcodes, false, lm, mol);
  }
  return 0;
}

int read_args(int argc, char const *argv[], layerModel &lm, collagenMolecule &mol)
{
  int errstate = 0;

  if(argc > 0)
  {
    for(int i = 0; i < argc; ++i)
      {
        int errcodes = 0;
        std::string strarg = argv[i];
        process_arg(strarg, &errcodes, true, lm, mol);

      }
  }
  return errstate;
}

int parse_all_args(int argc, char const *argv[], layerModel &lm, collagenMolecule &mol)
{
  std::cout << "\n# Parsing config file and command line arguments";
  std::string cpatharg = get_config_path(argc, argv);
  if (!cpatharg.empty()) {
    filePaths.configpath = cpatharg;
  }

  if (filePaths.configpath.empty()) {
    std::cout << "\n#  -> no config file selected, using defaults" << std::endl;
  }
  else {
    std::cout << "\n#  -> reading from config path " << filePaths.configpath << std::endl;
    std::cout << "#  -> Options:";
    read_config_file(filePaths.configpath, lm, mol);
  }

  int argerr = read_args(argc, argv, lm, mol);
  if (argerr > 0) {
    std::cout << "\n#  -> improper arguments supplied, exiting now" << std::endl;
    return 0;
  }
  else if (argerr < 0){
    return 0;
  }
  return 1;
}

void readAtomInfos(collagenMolecule &mol)
{
  if (flags.consoleOutput) {
    std::cout << "\n#\n#";
    std::cout << "\n# Reading atomic collagen information from file ";
    std::cout << filePaths.inputpath;
  }

  if(path_exists(filePaths.inputpath) == 1) {
    std::cout << "\n#  ERROR: path does not exist or could not be accessed: ";
    std::cout << filePaths.inputpath;
    std::cout << " -> using default";
    return;
  }

  int err = 0;
  std::string line;
  std::ifstream file;
  file.open(filePaths.inputpath, std::ios::in);
  std::getline(file, line);
  file.close();
  if (line.at(0) == '#') {
    line = line.erase(0, 1);
    line = remove_spaces(line);
    line = remove_formatting(line);
    if (line == "charges") {
      if (flags.consoleOutput) {
        std::cout << "\n#  -> file contains charge information";
      }
      mol.readCharges(filePaths.inputpath);
    } else if (line == "types") {
      if (flags.consoleOutput) {
        std::cout << "\n#  -> file contains type information";
      }
      mol.readTypes(filePaths.inputpath);
    } else {
      err = 2;
    }
  } else {
    err = 1;
  }

  if (err == 1) {
    std::cout << "\n# ERROR: Input file missing header file information!";
    return;
  }
  if (err == 2) {
    std::cout << "\n# ERROR: Input file has illegal header type '";
    std::cout << line << "'!";
    return;
  }

  if (flags.printAtomInfo && flags.consoleOutput) {
    mol.printAtoms();
  }
  if (flags.printMoleculeInfo && flags.consoleOutput) {
    mol.printMoleculeInfo();
  }
}

void prioParse(int argc, char const *argv[])
{
  for (int i = 0; i < argc; i++) {
    if (std::string(argv[i]) == "--s") {
      std::cout.setstate(std::ios_base::failbit);
    }
  }
}

void debugInfo(collagenMolecule mol)
{
  std::cout << "#";
  std::cout << "\n# Debug-Info:";


  std::cout << "\n#\t.Input filepaths:";

  std::cout << "\n#\t\t.molecule input: ";
  std::cout << filePaths.inputpath;


  std::cout << "\n#\t.Output filepaths:";

  std::cout << "\n#\t\t.general output: ";
  std::cout << filePaths.outputpath;

  std::cout << "\n#\t\t.md output folder: ";
  std::cout << filePaths.md_outputpath;

  std::cout << "\n#\t\t.md topology file: ";
  std::cout << mol.parametersMD.topFile;

  std::cout << "\n#\t\t.md lammps input file: ";
  std::cout << mol.parametersMD.inputFile;

  std::cout << "\n#\t\t.md lammps log filepath: ";
  std::cout << mol.parametersMD.logFile;

  std::cout << "\n#\t\t.md lammps dump filepath: ";
  std::cout << mol.parametersMD.dumpFile;
}

void programInfo()
{
  std::cout << "#";
  // std::cout << "\n# This is a program designed around the collagen type I";
  // std::cout << " molecule.";
  std::cout << "\n# run with option '-h' for help.";
  // std::cout << "\n# There are various functionalities.";
  // std::cout << "\n# See the readme, the config file, or run program with";
  // std::cout << " option \"-h\" for help.";
  // std::cout << "\n#";
  // std::cout << "\n# To use a config file, run program with option";
  // std::cout << " \"--config /path/to/config/file.config\".";
  // std::cout << "\n#";
  // std::cout << "\n#";
  // std::cout << "\n#";
  // std::cout << "\n# Running program...";
  std::cout << "\n#";
  std::cout << "\n#";
}

void printOptions(collagenMolecule mol, layerModel lm)
{
  if (flags.consoleOutput) {
    std::cout << "\n#\t.console output active";
  } else {
    std::cout << "\n#\t.console output inactive";
    return;
  }

  if (flags.measureTime) {
    std::cout << "\n#\t.measuring computation time";
  }

  if (flags.printMoleculeInfo) {
    std::cout << "\n#\t.printing molecule information";
  }
  if (flags.printAtomInfo) {
    std::cout << "\n#\t.printing atomic information";
  }

  if (flags.layermodel) {
    std::cout << "\n#\t.original output active";
    if (lm.parameters.chargehash) {
      std::cout << "\n#\t  option hashed output";
    }
    if (lm.parameters.xyz_outputs) {
      std::cout << "\n#\t  option xyz output";
    }
    if (filePaths.outputpath.find(filePaths.file_extension) != std::string::npos)
    {
        filePaths.outputpath = filePaths.outputpath.substr(0, filePaths.outputpath.find(filePaths.file_extension));
    }
    if (lm.parameters.csv_output) {
      filePaths.outputpath += filePaths.csv_extension;
    } else {
      filePaths.outputpath += filePaths.file_extension;
    }
    std::cout << "\n#\t  writing to: " << filePaths.outputpath;
  }

  if (flags.md) {
    std::cout << "\n#\t.MD active";
    if (mol.parametersMD.script) {
      std::cout << "\n#\t  option bash scripts";
    }
    if (mol.parametersMD.topology) {
      std::cout << "\n#\t  option topology";
    }
    if (mol.parametersMD.input) {
      std::cout << "\n#\t  option LAMMPS input";
    }
  }

  if (flags.development) {
    std::cout << "\n#\t.running developmental part";
  }
}

void getTime(time_point &tp)
{
  tp = std::chrono::high_resolution_clock::now();
}

void printCompTime(time_point start, time_point end)
{
  std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  if (flags.measureTime && flags.consoleOutput) {
    std::cout << "\n#\n#\n# Run time: " << duration.count() << " s.";
  }
}

void print_help()
{
  std::cout << "\n\nHelp to be done!";
}
