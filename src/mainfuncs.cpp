#include "mainfuncs.hpp"


/* Variables */
overrides_ overrides;
extern filePaths_ filePaths;
extern flags_ flags;

/* Functions */
int process_arg(std::string strarg, int* errcodes, bool dashes, collagenFibril &fib)
{
  //io overrides
    int errstate = 0;

    std::string d = dashes ? "-" : "";
    std::string dd = dashes ? "--" : "";
    if(overrides.input_override)
    {
      if(verify_path(strarg,errcodes) == 0) {
        filePaths.inputpath = strarg;
      }
      overrides.input_override = false;
    }

    if(overrides.charges_input_override)
    {
      if(verify_path(strarg,errcodes) == 0) {
        filePaths.charges_inputpath = strarg;
      }
      overrides.charges_input_override = false;
    }

    if(overrides.types_input_override)
    {
      if(verify_path(strarg,errcodes) == 0) {
        filePaths.types_inputpath = strarg;
      }
      overrides.types_input_override = false;
    }

    if(overrides.output_override)
    {
      filePaths.outputpath = strarg;
      overrides.output_override = false;
    }

    if(overrides.md_output_override)
    {
      filePaths.md_outputpath = strarg;
      overrides.md_output_override = false;
    }

    //spatial overrides

    if(overrides.layers_override)
    {
      fib.layers = safe_read_integer(fib.layers,strarg,errcodes);
      overrides.layers_override = false;
    }
    if(overrides.diameter_override)
    {
      fib.mol.diameterAtom = safe_read_double(fib.mol.diameterAtom,strarg,errcodes);
      overrides.diameter_override = false;
    }
    if(overrides.distance_override)
    {
      fib.mol.distanceAtoms = safe_read_double(fib.mol.distanceAtoms,strarg,errcodes);
      overrides.distance_override = false;
    }

    //potential overrides

    if(overrides.ljsteps_override)
    {
      fib.parameters.lj_steps = safe_read_integer(fib.parameters.lj_steps,strarg,errcodes);
      overrides.ljsteps_override = false;
    }
    if(overrides.ljstepsize_override)
    {
      fib.parameters.lj_stepsize = safe_read_double(fib.parameters.lj_stepsize,strarg,errcodes);
      overrides.ljstepsize_override = false;
    }
    if(overrides.cdsteps_override)
    {
      fib.parameters.cd_steps = safe_read_integer(fib.parameters.cd_steps,strarg,errcodes);
      overrides.cdsteps_override = false;
    }
    if(overrides.cdstepsize_override)
    {
      fib.parameters.cd_stepsize = safe_read_double(fib.parameters.cd_stepsize,strarg,errcodes);
      overrides.cdstepsize_override = false;
    }
    if(overrides.ljmin_override)
    {
      fib.parameters.lj_min = safe_read_double(fib.parameters.lj_min,strarg,errcodes);
      overrides.ljmin_override = false;
    }
    if(overrides.cdmin_override)
    {
      fib.parameters.cd_min = safe_read_double(fib.parameters.cd_min,strarg,errcodes);
      overrides.cdmin_override = false;
    }
    if(overrides.ljcut_override)
    {
      fib.parameters.lj_cutoff = safe_read_double(fib.parameters.lj_cutoff,strarg,errcodes);
      overrides.ljcut_override = false;
    }
    if(overrides.cdcut_override)
    {
      fib.parameters.cd_cutoff = safe_read_double(fib.parameters.cd_cutoff,strarg,errcodes);
      overrides.cdcut_override = false;
    }

    // MD related
    if(overrides.MD_walltime_override)
    {
      fib.parametersMD.walltime = safe_read_integer(fib.parametersMD.walltime,strarg,errcodes);;
      overrides.MD_walltime_override = false;
    }
    if(overrides.MD_cores_override)
    {
      fib.parametersMD.cores = safe_read_integer(fib.parametersMD.cores,strarg,errcodes);;
      overrides.MD_cores_override = false;
    }
    if(overrides.MD_numMolperDim_override)
    {
      fib.parametersMD.numMolperDim = safe_read_integer(fib.parametersMD.numMolperDim,strarg,errcodes);;
      overrides.MD_numMolperDim_override = false;
    }
    if(overrides.MD_dielectric_override)
    {
      fib.parametersMD.dielectric = safe_read_double(fib.parametersMD.dielectric,strarg,errcodes);
      overrides.MD_dielectric_override = false;
    }
    if(overrides.MD_LJepsilon_override)
    {
      fib.parametersMD.LJepsilon = safe_read_double(fib.parametersMD.LJepsilon,strarg,errcodes);
      overrides.MD_LJepsilon_override = false;
    }
    if(overrides.MD_LJcutoff_override)
    {
      fib.parametersMD.lj_cutoff = safe_read_double(fib.parametersMD.lj_cutoff,strarg,errcodes);
      overrides.MD_LJcutoff_override = false;
    }
    if(overrides.MD_CDcutoff_override)
    {
      fib.parametersMD.cd_cutoff = safe_read_double(fib.parametersMD.cd_cutoff,strarg,errcodes);
      overrides.MD_CDcutoff_override = false;
    }
    if(overrides.MD_timestep_override)
    {
      fib.parametersMD.timestep = safe_read_double(fib.parametersMD.timestep,strarg,errcodes);
      overrides.MD_timestep_override = false;
    }
    if(overrides.MD_runtime_override)
    {
      fib.parametersMD.runtime = safe_read_integer(fib.parametersMD.runtime,strarg,errcodes);
      overrides.MD_runtime_override = false;
    }

    //error handling

    if (parse_errs(*errcodes,strarg,true) != 0)
    {
      errstate = 1;
    }

    //bool flags and help

    if(flag(strarg,d+"h") || flag(strarg,dd+"help"))
    {
      print_help();
      return -1;
    }

    if(flag(strarg,d+"c") || flag(strarg,dd+"chargehash"))
    {
      flags.charge_hashed_outputs = true;
    }

    if(flag(strarg,d+"x") || flag(strarg,dd+"xyz"))
    {
      flags.xyz_outputs = true;
    }

    if(flag(strarg,d+"rc") || flag(strarg,dd+"readCharges"))
    {
      flags.readCharges = true;
    }

    if(flag(strarg,d+"rt") || flag(strarg,dd+"readTypes"))
    {
      flags.readTypes = true;
    }

    if(flag(strarg,d+"pai") || flag(strarg,dd+"printAtomInfo"))
    {
      flags.printAtomInfo = true;
    }

    if(flag(strarg,d+"moli") || flag(strarg,dd+"printMoleculeInfo"))
    {
      flags.printMoleculeInfo = true;
    }

    if(flag(strarg,d+"anneo") || flag(strarg,dd+"annesOutput"))
    {
      flags.annesOutput = true;
    }

    if(flag(strarg,d+"time") || flag(strarg,dd+"measureTime"))
    {
      flags.measureTime = true;
    }

    if(flag(strarg,d+"dev") || flag(strarg,dd+"development"))
    {
      flags.development = true;
    }


    //io overrides

    if(flag(strarg,d+"i") || flag(strarg,dd+"input"))
    {
      overrides.input_override = true;
    }
    if(flag(strarg,d+"ci") || flag(strarg,dd+"chargesinput"))
    {
      overrides.charges_input_override = true;
    }
    if(flag(strarg,d+"ti") || flag(strarg,dd+"typesinput"))
    {
      overrides.types_input_override = true;
    }

    if(flag(strarg,d+"o") || flag(strarg,dd+"output"))
    {
      overrides.output_override = true;
    }
    if(flag(strarg,d+"mdo") || flag(strarg,dd+"mdoutput"))
    {
      overrides.md_output_override = true;
    }

    //spatial overrides

    if(flag(strarg,d+"l") || flag(strarg,dd+"layers"))
    {
      overrides.layers_override = true;
    }
    if(flag(strarg,d+"dia") || flag(strarg,dd+"diameter"))
    {
      overrides.diameter_override = true;
    }
    if(flag(strarg,d+"d") || flag(strarg,dd+"dist"))
    {
      overrides.distance_override = true;
    }

    //potential overrides

    if(flag(strarg,d+"slj") || flag(strarg,dd+"ljsteps"))
    {
      overrides.ljsteps_override = true;
    }
    if(flag(strarg,d+"sslj") || flag(strarg,dd+"ljstepssize"))
    {
      overrides.ljstepsize_override = true;
    }
    if(flag(strarg,d+"scd") || flag(strarg,dd+"cdsteps"))
    {
      overrides.cdsteps_override = true;
    }
    if(flag(strarg,d+"sscd") || flag(strarg,dd+"cdstepssize"))
    {
      overrides.cdstepsize_override = true;
    }
    if(flag(strarg,d+"mcd") || flag(strarg,dd+"cdmin"))
    {
      overrides.cdmin_override = true;
    }
    if(flag(strarg,d+"mlj") || flag(strarg,dd+"ljmin"))
    {
      overrides.ljmin_override = true;
    }


    if(flag(strarg,d+"ccd") || flag(strarg,dd+"cdcut"))
    {
      overrides.cdcut_override = true;
    }
    if(flag(strarg,d+"clj") || flag(strarg,dd+"ljcut"))
    {
      overrides.ljcut_override = true;
    }

    // MD related
    if(flag(strarg,d+"md_s") || flag(strarg,dd+"MD_script"))
    {
      fib.parametersMD.outputScript = true;
    }
    if(flag(strarg,d+"md_t") || flag(strarg,dd+"MD_topology"))
    {
      fib.parametersMD.outputTopology = true;
    }
    if(flag(strarg,d+"md_li") || flag(strarg,dd+"MD_LAMMPS_input"))
    {
      fib.parametersMD.outputLAMMPSinput = true;
    }
    if(flag(strarg,d+"md_sb"))
    {
      fib.parametersMD.scriptbuild = true;
    }
    if(flag(strarg,d+"md_rigid") || flag(strarg,dd+"MD_rigid"))
    {
      fib.parametersMD.rigid = true;
    }
    if(flag(strarg,d+"md_wt") || flag(strarg,dd+"MD_walltime"))
    {
      overrides.MD_walltime_override = true;
    }
    if(flag(strarg,d+"md_c") || flag(strarg,dd+"MD_cores"))
    {
      overrides.MD_cores_override = true;
    }
    if(flag(strarg,d+"md_N") || flag(strarg,dd+"MD_N"))
    {
      overrides.MD_numMolperDim_override = true;
    }
    if(flag(strarg,d+"md_die") || flag(strarg,dd+"MD_dielectric"))
    {
      overrides.MD_dielectric_override = true;
    }
    if(flag(strarg,d+"md_lje") || flag(strarg,dd+"MD_LJepsilon"))
    {
      overrides.MD_LJepsilon_override = true;
    }
    if(flag(strarg,d+"md_ljc") || flag(strarg,dd+"MD_LJcutoff"))
    {
      overrides.MD_LJcutoff_override = true;
    }
    if(flag(strarg,d+"md_cdc") || flag(strarg,dd+"MD_CDcutoff"))
    {
      overrides.MD_CDcutoff_override = true;
    }
    if(flag(strarg,d+"md_ts") || flag(strarg,dd+"MD_timestep"))
    {
      overrides.MD_timestep_override = true;
    }
    if(flag(strarg,d+"md_rt") || flag(strarg,dd+"MD_runtime"))
    {
      overrides.MD_runtime_override = true;
    }

    return errstate;
}

int read_config_file(std::string path, collagenFibril &fib)
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
          process_arg(token, &errcodes, false, fib);
          s.erase(0, pos + delimiter.length());
      }
      process_arg(s, &errcodes, false, fib);
  }
  return 0;
}

int read_args(int argc, char const *argv[], collagenFibril &fib)
{
  int errstate = 0;

  if(argc > 0)
  {
    for(int i = 0; i < argc; ++i)
      {
        int errcodes = 0;
        std::string strarg = argv[i];
        process_arg(strarg, &errcodes, true, fib);

      }
  }
  return errstate;
}

int parse_all_args(int argc, char const *argv[], collagenFibril &fib)
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
    std::cout << "#  -> Options:" << std::endl;
    read_config_file(filePaths.configpath, fib);
  }

  int argerr = read_args(argc, argv, fib);
  if (argerr > 0) {
    std::cout << "\n#  -> improper arguments supplied, exiting now" << std::endl;
    return 0;
  }
  else if (argerr < 0){
    return 0;
  }
  return 1;
}

void readAtomInfos(collagenFibril &fib)
{
  if (flags.readCharges) {
    fib.mol.readCharges(filePaths.charges_inputpath);
  }
  if (flags.readTypes) {
    fib.mol.readTypes(filePaths.types_inputpath);
  }
  if (!flags.readCharges && !flags.readTypes) {
    std::cout << "#\n  WARNING: No collagen molecule information has been read!";
    return;
  } else if (flags.printAtomInfo) {
    fib.mol.printAtoms();
  }
  if (flags.printMoleculeInfo) {
    fib.printMoleculeInfo();
  }
}

void programInfo()
{
  std::cout << "#";
  std::cout << "\n# This is a program designed around the collagen type I";
  std::cout << " molecule.";
  std::cout << "\n# There are various functionalities.";
  std::cout << "\n# See the readme, the config file, or run program with";
  std::cout << " option \"-h\" for help.";
  std::cout << "\n#";
  std::cout << "\n# To use a config file, run program with option";
  std::cout << " \"--config /path/to/config/file.config\".";
  std::cout << "\n#";
  std::cout << "\n#";
  std::cout << "\n#";
  std::cout << "\n# Running program...";
  std::cout << "\n#";
  std::cout << "\n#";
  std::cout << "\n#";
}

void printOptions(collagenFibril fib)
{
  if (flags.readTypes) {
    std::cout << "\n#\treading atom types from " << filePaths.types_inputpath;
  } else if (flags.readCharges) {
    std::cout << "\n#\treading atom charges from " << filePaths.charges_inputpath;
  }
  if (flags.printMoleculeInfo) {
    std::cout << "\n#\tprinting molecule information";
  }
  if (flags.printAtomInfo) {
    std::cout << "\n#\tprinting atomic information";
  }
  if (flags.annesOutput) {
    std::cout << "\n#\tcreating Anne's output";
    if (flags.charge_hashed_outputs) {
      std::cout << "\n#\t  option hashed output";
    }
    if (flags.xyz_outputs) {
      std::cout << "\n#\t  option xyz output";
    }
    std::cout << "\n#\t  writing to: " << filePaths.outputpath;
  }
  if (flags.development) {
    std::cout << "\n#\trunning developmental part";
    if (fib.parametersMD.outputTopology) {
      std::cout << "\n#\t  generating LAAMPS topology file for ";
      std::cout << pow(fib.parametersMD.numMolperDim, 3) << " molecules";
      std::cout << " in " << filePaths.md_outputpath;
    }
    if (fib.parametersMD.outputLAMMPSinput) {
      std::cout << "\n#\t  generating input file for LAMMPS";
      std::cout << " in " << filePaths.md_outputpath;
    }
    if (fib.parametersMD.outputScript) {
      std::cout << "\n#\t  generating bash scripts for LAMMPS content";
      std::cout << " in " << filePaths.md_outputpath;
    }
  }
  std::cout << "\n#";
}

void print_help()
{
  std::cout << "\n\nHelp to be done!";
}
