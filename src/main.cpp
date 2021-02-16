#include "main.hpp"
#include "mainfuncs.hpp"
#include "parse.hpp"
#include "collmol.hpp"
#include "collfibril.hpp"
#include "random.hpp"

using namespace std;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point;


/* Settings */
filePaths_ filePaths;       /* IO File paths */
flags_ flags;               /* Various flags, see main.hpp */

int main(int argc, char const *argv[])
{

  collagenFibril fib;

  if (parse_all_args(argc, argv, fib) != 1) {
    return 0;
  }

  time_point start = chrono::high_resolution_clock::now();
  cout << ":: Computation start";

  /************************************************************************/

  cout << "\n\n:: running computation...";
  /* Only global Emin */
  fib.mol.readAtoms(filePaths.inputpath);
  fib.singleEmin();

  /************************************************************************/

  time_point end = chrono::high_resolution_clock::now();
  chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
  cout << "\n\n\n:: ... finished in " << duration.count() << " s.\n\n\n";

  /************************************************************************/

  return 0;
}
