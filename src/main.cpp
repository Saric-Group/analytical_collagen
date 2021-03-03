#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"
#include "random.hpp"
#include "xyz.hpp"
#include "funcs.hpp"
#include "md.hpp"

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;


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

  time_point start;
  if (flags.measureTime) {
    std::cout << "#\tmeasuring computation time";
    start = std::chrono::high_resolution_clock::now();
  }
  printOptions();

  readAtomInfos(fib);

  if (flags.annesOutput) {
    fib.singleEmin();
  }

  /************************************************************************/

  if (flags.measureTime) {
    time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "\n#\n#\n# Run time: " << duration.count() << " s.";
  }

  /************************************************************************/

  std::cout << "\n#\n# Exiting...\n#\n";
  return 0;

}
