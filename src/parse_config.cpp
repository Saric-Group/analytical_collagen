#include "parse_config.hpp"


/* Functions */
void parse_config()
{
  setlocale(LC_ALL, "");

  const char *remoteConfig = "exec#curl -sS https://github.com/Saric-Group/analytical_collagen/tree/new_config%2Bcleanup/config/test.config";

  const char *scope = "collagen";

  bool consOut;

  config4cpp::Configuration *cfg = config4cpp::Configuration::create();

  try {
    cfg->parse(remoteConfig);
    consOut = cfg->lookupBoolean(scope, "consoleOutput");
  } catch (const config4cpp::ConfigurationException &ex) {
    std::cerr << ex.c_str() << "\n";
    cfg->destroy();
  }
  std::cout << "\n# -> console output: " << consOut;
  cfg->destroy();
}
