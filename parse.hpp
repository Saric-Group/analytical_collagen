#include "coll.h"

using namespace std;

//error codes

int misc_err = 1;
int d_parse_err = 2;
int d_limit_err = 4;
int i_parse_err = 8;
int i_limit_err = 16;
int path_exist_err = 32;

double read_double(string arg, int* err){
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

int read_integer(string arg, int* err){
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

int parse_errs(const int err, const string arg, const bool verbose){
  int errstate = 0;
  if(err&d_parse_err){
    errstate = 1;
    if(verbose){
      cout << "error parsing argument: " << arg << " as double" << endl;
    }
  }
  if(err&d_limit_err){
    errstate = 1;
    if(verbose){
      cout << "argument overflowed buffer: " << arg << " as double" << endl;
    }
  }
  if(err&i_parse_err){
    errstate = 1;
    if(verbose){
      cout << "error parsing argument: " << arg << " as integer" << endl;
    }
  }
  if(err&i_limit_err){
    errstate = 1;
    if(verbose){
      cout << "argument overflowed buffer: " << arg << " as integer" << endl;
    }
  }
  if(err&path_exist_err){
    errstate = 1;
    if(verbose){
      cout << "path does not exist or could not be accessed: " << arg << endl;
    }
  }
  return errstate;
}

int path_exists(const std::string& name){
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return 0;
    } else {
        return 1;
    }   
}

int verify_path(const string path, int* err){
  if(path_exists(path)!=0){
    *err += path_exist_err;
  }
  return *err >0 ? 1 : 0;
}

string get_config_path(const int argc, char const *argv[])
{
  int errcodes = 0;
  int configerr = -1; // -1 for no config file, 0 for succesful read, 1 for errors
  bool has_config = false;
  string cpath = "";
  if(argc > 0)
  {                                                                          
    for(int i=0; i<argc; ++i)
      {
        string strarg = argv[i];
        if(has_config){
          configerr = 0;
          cpath = strarg;
          has_config = false;
        }
        if(strarg.compare("--conf")==0 || strarg.compare("--config")==0 || strarg.compare("-f")==0){
          has_config = true;
        }
      }
  }
  if(cpath.compare("") != 0){
    int verifyerr = verify_path(cpath,&errcodes);
    configerr = parse_errs(errcodes,cpath,true);
  }

  return configerr==0?cpath:"";
}

bool flag(string arg,string str){
  return arg.compare(str)==0;
}

string remove_spaces(string input)
{
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
  return input;
}

int safe_read_integer(int orig,string strarg,int* errcodes){
  int parsed = read_integer(strarg,errcodes);
  return *errcodes == 0 ? parsed : orig;
}

double safe_read_double(double orig,string strarg,int* errcodes){
  double parsed = read_double(strarg,errcodes);
  return *errcodes == 0 ? parsed : orig;
}