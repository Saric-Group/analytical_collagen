#ifndef PARSE_CONFIG
#define PARSE_CONFIG

#include "main.hpp"
#include "collmol.hpp"
#include "layermodel.hpp"
#include <locale.h>
#include "config4cpp/include/config4cpp/Configuration.h"
#include "config4cpp/include/config4cpp/SchemaValidator.h"
// #include "FallbackConfig.h"


/* Functions */
int parse(int argc, char const *argv[], collagenMolecule &mol, layerModel &lm);

#endif
