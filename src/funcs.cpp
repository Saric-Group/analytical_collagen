#include "funcs.hpp"


/* Functions */
void progresssBar(double progress, std::string msg)
{
  int barwidth = 50;
  int pos = 1.0 * barwidth * progress;
  std::cout << "\r# " << msg << " [";
  for (int i = 0; i < barwidth; i++) {
    if (i < pos) {
      std::cout << "#";
    } else {
      std::cout << "-";
    }
  }
  std::cout << "] ";
  std::cout << std::setw(5);
  std::cout << std::fixed;
  std::cout << std::setprecision(1);
  std::cout << 100.0 * progress << "%";
  std::cout.flush();
  std::cout << std::setprecision(4);
}
