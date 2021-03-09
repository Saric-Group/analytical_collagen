#ifndef PARSE
#define PARSE

#include "main.hpp"


/* Functions */
std::string lastchar(std::string input);
std::string clean_to_string(double f);
std::string replace_char(std::string str, char ch1, char ch2);
double read_double(std::string arg, int* err);
int read_integer(std::string arg, int* err);
int parse_errs(const int err, const std::string arg, const bool verbose);
int path_exists(const std::string& name);
int verify_path(const std::string path, int* err);
std::string get_config_path(const int argc, char const *argv[]);
bool flag(std::string arg, std::string str);
std::string remove_spaces(std::string input);
std::string remove_formatting(std::string input);
int safe_read_integer(int orig, std::string strarg, int* errcodes);
double safe_read_double(double orig, std::string strarg, int* errcodes);

#endif
