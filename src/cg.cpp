#include "cg.hpp"


/* Functions */
std::vector<double> smoothenSMCA(int delta,
                                 std::vector<double> vec,
                                 bool avg /*= false*/)
{
  double sma = 0;
  double help;
  std::vector<double> svec;

  for (int i = 0; i <= delta; i++) {
    help = vec[i];
    if (avg) help /= 2 * delta + 1;
    sma += help;
  }
  svec.push_back(sma);

  for (int i = 1; i < (int) vec.size(); i++) {
    if (i - 1 - delta < 0) {
      help = vec[i + delta];
      if (avg) help /= 2 * delta + 1;
      sma += help;
    } else if (i + delta > (int) vec.size() - 1) {
      help = vec[i - 1 - delta];
      if (avg) help /= 2 * delta + 1;
      sma -= help;
    } else {
      help = vec[i + delta] - vec[i - 1 - delta];
      if (avg) help /= 2 * delta + 1;
      sma += help;
    }
    svec.push_back(sma);
  }

  return svec;
}

std::vector<double> binNormalize(int n, std::vector<double> vec)
{
  std::vector<double> nvec;
  double ma;
  int binSize = (int) floor(vec.size() / n);
  int split = vec.size() % n;
  for (int i = 0; i < split; i++) {
    // cout << "\nBin " << i + 1 << " adding:\n";
    ma = 0;
    for (int j = 0; j <= binSize; j++) {
      // cout << i * (binSize + 1) + j << " ";
      ma += vec[i * (binSize + 1) + j];
    }
    nvec.push_back(ma);
  }
  for (int i = split; i < n - 1; i++) {
    // cout << "\nBin " << i + 1 << " adding:\n";
    ma = 0;
    for (int j = 0; j < binSize; j++) {
      // cout << i * binSize + split + j << " ";
      ma += vec[i * binSize + split + j];
    }
    nvec.push_back(ma);
  }
  ma = 0;
  // cout << "\nBin " << n << " adding:\n";
  for (int i = (n - 1) * binSize + split; i < (int) vec.size(); i++) {
    // cout << i << " ";
    ma += vec[i];
  }
  nvec.push_back(ma);

  // cout << "\n\ntest: " << vec.size() % n;

  double sum = 0.0;
  for (int i = 0; i < (int) nvec.size(); i++) {
    sum += nvec[i];
  }
  // cout << "\nSize: " << nvec.size();
  // cout << "\nIntegral: " << sum;

  return nvec;
}
