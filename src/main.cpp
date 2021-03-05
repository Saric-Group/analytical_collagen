#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"
#include "random.hpp"
#include "xyz.hpp"
#include "funcs.hpp"
#include "md.hpp"
#include "cg.hpp"


/* Settings */
filePaths_ filePaths;       /* IO File paths */
flags_ flags;               /* Various flags, see main.hpp */


int main(int argc, char const *argv[])
{
  programInfo();

  collagenFibril fib;

  if (parse_all_args(argc, argv, fib) != 1) {
    return 0;
  }

  time_point start, end;
  getTime(start);

  printOptions(fib);

  readAtomInfos(fib);

  if (flags.originalOutput) {
    fib.singleEmin();
  }

  if (flags.mdOutput) {
    createLAMMPSfiles(fib);
  }

  if (flags.development) {
  }

  /************************************************************************/

  getTime(end);
  printCompTime(start, end);

  /************************************************************************/

  std::cout << "\n#\n# exiting...\n#\n";

  return 0;

}
