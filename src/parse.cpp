#include "parse.hpp"


/* Variables */
//error codes
int misc_err = 1;
int d_parse_err = 2;
int d_limit_err = 4;
int i_parse_err = 8;
int i_limit_err = 16;
int path_exist_err = 32;


/* Functions */
std::string lastchar(std::string input)
{
  return input.substr(input.size() - 1);
}

std::string clean_to_string(double f)
{
  std::string str = std::to_string (f);
  str.erase(str.find_last_not_of('0') + 1, std::string::npos);
  if(lastchar(str) == ".") {
    str.erase(str.find_last_not_of('.') + 1, std::string::npos);
  }
  return str;
}

std::string replace_char(std::string str, char ch1, char ch2)
{
  for (int i = 0; i < (int) str.length(); ++i) {
    if (str[i] == ch1)
      str[i] = ch2;
  }
  return str;
}

double read_double(std::string arg, int* err)
{
  double argd = 0.0;
  try {
    argd = std::stod(arg);
  } catch (const std::invalid_argument &) {
    *err += d_parse_err;
  } catch (const std::out_of_range &) {
    *err += d_limit_err;
  }
  return argd;
}

int read_integer(std::string arg, int* err)
{
  int argi = 0;
  try {
    argi = std::stoi(arg);
  } catch (const std::invalid_argument &) {
    *err += i_parse_err;
  } catch (const std::out_of_range &) {
    *err += i_limit_err;
  }
  return argi;
}

int parse_errs(const int err, const std::string arg, const bool verbose)
{
  int errstate = 0;
  if(err&d_parse_err){
    errstate = 1;
    if(verbose){
      std::cout << "#  ERROR: parsing argument: " << arg << " as double" << std::endl;
    }
  }
  if(err&d_limit_err){
    errstate = 1;
    if(verbose){
      std::cout << "#  ERROR: argument overflowed buffer: " << arg << " as double" << std::endl;
    }
  }
  if(err&i_parse_err){
    errstate = 1;
    if(verbose){
      std::cout << "#  ERROR: error parsing argument: " << arg << " as integer" << std::endl;
    }
  }
  if(err&i_limit_err){
    errstate = 1;
    if(verbose){
      std::cout << "#  ERROR: argument overflowed buffer: " << arg << " as integer" << std::endl;
    }
  }
  if(err&path_exist_err){
    errstate = 1;
    if(verbose){
      std::cout << "\n#  ERROR: path does not exist or could not be accessed: " << arg << std::endl;
    }
  }
  return errstate;
}

int path_exists(const std::string& name)
{
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return 0;
    } else {
        return 1;
    }
}

int verify_path(const std::string path, int* err)
{
  if(path_exists(path) != 0) {
    *err += path_exist_err;
  }
  return *err > 0 ? 1 : 0;
}

std::string get_config_path(const int argc, char const *argv[])
{
  int errcodes = 0;
  int configerr = -1; // -1 for no config file, 0 for succesful read, 1 for errors
  bool has_config = false;
  std::string cpath = "";
  if(argc > 0)
  {
    for(int i = 0; i < argc; ++i)
      {
        std::string strarg = argv[i];
        if(has_config) {
          configerr = 0;
          cpath = strarg;
          has_config = false;
        }
        if(strarg.compare("--conf") == 0 || strarg.compare("--config") == 0 || strarg.compare("-f") == 0) {
          has_config = true;
        }
      }
  }
  if(cpath.compare("") != 0) {
    verify_path(cpath,&errcodes);
    configerr = parse_errs(errcodes,cpath,true);
  }

  return configerr == 0 ? cpath : "";
}

bool flag(std::string arg, std::string str)
{
  return arg.compare(str) == 0;
}

std::string remove_spaces(std::string input)
{
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
  return input;
}

int safe_read_integer(int orig, std::string strarg, int* errcodes)
{
  int parsed = read_integer(strarg,errcodes);
  return *errcodes == 0 ? parsed : orig;
}

double safe_read_double(double orig, std::string strarg, int* errcodes)
{
  double parsed = read_double(strarg, errcodes);
  return *errcodes == 0 ? parsed : orig;
}
