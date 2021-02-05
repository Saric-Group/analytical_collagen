#include "mainfuncs.hpp"


/* Variables */
overrides_ overrides;
extern parameters_ parameters;
extern filePaths_ filePaths;
extern flags_ flags;

/* Functions */
std::string hash_output(collagenFibril &fib)
{
  //create a unique code for the output file

    unsigned long f_charges = 0;

    std::string fflags = "";

    unsigned long lulen = sizeof(unsigned long)*8;

    int fcount = 0;

    double chargesum=0.0;
    int chargecount = 0;

    for(int c = 0; c < (int) fib.mol.charges.size(); c++)
    {
      if(fib.mol.charges[c] != 0){
        chargesum += fib.mol.charges[c];
        chargecount++;
      }
    }

    double chargemean = chargesum / ((double) chargecount);

    fflags += "-";
    fflags += std::to_string(fib.mol.numAtoms);
    fflags += "-";
    fflags += std::to_string(fib.mol.distanceAtoms);
    fflags += "-";
    fflags += std::to_string(chargesum);
    fflags += "-";
    fflags += std::to_string(chargemean);
    fflags += "-";



    for(int i = 0; i < fib.mol.numAtoms; i++)
    {
      if(fcount >= lulen){
        fflags += std::to_string(f_charges);
        fflags += "-";
        f_charges = 0;
        fcount = 0;
      }
      if((int) fib.mol.charges[i] != 0)
      {
        f_charges += pow(2,(int) (lulen - fcount - 1));
      }

      fcount++;
    }

    fflags += std::to_string(f_charges);


    if (filePaths.outputpath.find(filePaths.file_extension) != std::string::npos)
    {
      filePaths.outputpath = filePaths.outputpath.substr(0, filePaths.outputpath.find(filePaths.file_extension));
    }

    filePaths.outputpath += fflags;
    return filePaths.outputpath;
}

int process_arg(std::string strarg, int* errcodes, bool dashes, collagenFibril &fib)
{
  //io overrides
    int errstate = 0;

    std::string d = dashes ? "-" : "";
    std::string dd = dashes ? "--" : "";
    if(overrides.input_override)
    {
      if(verify_path(strarg,errcodes)==0){
        filePaths.inputpath = strarg;
      }
      overrides.input_override = false;
    }
    if(overrides.output_override)
    {
      filePaths.outputpath = strarg;
      overrides.output_override = false;
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
      parameters.lj_steps = safe_read_integer(parameters.lj_steps,strarg,errcodes);
      overrides.ljsteps_override = false;
    }
    if(overrides.ljstepsize_override)
    {
      parameters.lj_stepsize = safe_read_double(parameters.lj_stepsize,strarg,errcodes);
      overrides.ljstepsize_override = false;
    }
    if(overrides.cdsteps_override)
    {
      parameters.cd_steps = safe_read_integer(parameters.cd_steps,strarg,errcodes);
      overrides.cdsteps_override = false;
    }
    if(overrides.cdstepsize_override)
    {
      parameters.cd_stepsize = safe_read_double(parameters.cd_stepsize,strarg,errcodes);
      overrides.cdstepsize_override = false;
    }
    if(overrides.ljmin_override)
    {
      parameters.lj_min = safe_read_double(parameters.lj_min,strarg,errcodes);
      overrides.ljmin_override = false;
    }
    if(overrides.cdmin_override)
    {
      parameters.cd_min = safe_read_double(parameters.cd_min,strarg,errcodes);
      overrides.cdmin_override = false;
    }
    if(overrides.ljcut_override)
    {
      parameters.lj_cutoff = safe_read_double(parameters.lj_cutoff,strarg,errcodes);
      overrides.ljcut_override = false;
    }
    if(overrides.cdcut_override)
    {
      parameters.cd_cutoff = safe_read_double(parameters.cd_cutoff,strarg,errcodes);
      overrides.cdcut_override = false;
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


    //io overrides

    if(flag(strarg,d+"i") || flag(strarg,dd+"input"))
    {
      overrides.input_override = true;
    }
    if(flag(strarg,d+"o") || flag(strarg,dd+"output"))
    {
      overrides.output_override = true;
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


    if(flag(strarg,"-ccd") || flag(strarg,dd+"cdcut"))
    {
      overrides.cdcut_override = true;
    }
    if(flag(strarg,"-clj") || flag(strarg,dd+"ljcut"))
    {
      overrides.ljcut_override = true;
    }
    return errstate;
}

int read_config_file(std::string path, collagenFibril &fib)
{
  int errcodes = 0;
  std::ifstream infile("default.config");
  std::string line;
  while (std::getline(infile, line))
  {
      errcodes = 0;
      std::string s = remove_spaces(line);
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
  std::string cpatharg = get_config_path(argc, argv);
  if (!cpatharg.empty()) {
    filePaths.configpath = cpatharg;
  }

  if (filePaths.configpath.empty()) {
    std::cout << "no config file selected, using defaults" << std::endl;
  }
  else {
    std::cout << "reading from config path " << filePaths.configpath << std::endl;
    read_config_file(filePaths.configpath, fib);
  }

  int argerr = read_args(argc, argv, fib);
  if (argerr > 0) {
    std::cout << "improper arguments supplied, exiting now" << std::endl;
    return 0;
  }
  else if (argerr < 0){
    return 0;
  }
  return 1;
}
