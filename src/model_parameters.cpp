#include "model_parameters.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "pll-utils.hpp"

namespace pt { namespace pll {

ModelParameters ParseRaxmlInfo(const std::string& path)
{
  std::ifstream file(path);
  std::string read;
  std::string contents;

  while (std::getline(file, read)) {
    contents += read;
    contents.push_back('\n');
  }

  // initialize the array of base frequencies
  std::size_t pos1 = contents.find("frequencies: ");
  std::size_t pos2 = contents.find('\n', pos1);
  std::string sstr = contents.substr(pos1 + 13, pos2 - pos1 - 13);

  std::vector<std::string> freqvector = ssplit(sstr, ' ');
  std::vector<double> frequencies(freqvector.size());

  for (unsigned int i = 0; i < freqvector.size(); i++)
    frequencies[i] = std::stod(freqvector.at(i));

  // initialize the array of subst_params
  pos1 = contents.find("ac ag at cg ct gt:");
  pos2 = contents.find('\n', pos1);
  sstr = contents.substr(pos1 + 19, pos2 - pos1 - 19);
  std::vector<std::string> ratevector = ssplit(sstr, ' ');
  std::vector<double> subst_params(ratevector.size());

  for (unsigned int i = 0; i < ratevector.size(); i++)
    subst_params[i] = std::stod(ratevector.at(i));
  if (ratevector.size() != (((freqvector.size()) * freqvector.size() - 4) / 2)) {
    throw std::invalid_argument("Wrong number of rate values.");
  }

  // we'll use RATE_CATS rate categories, and currently initialize them to 0
  std::vector<double> rate_cats(RATE_CATS, 0.0);

  // initialize the alpha value
  pos1 = contents.find("alpha[0]: ");
  pos2 = contents.find(' ', pos1 + 10);
  sstr = contents.substr(pos1 + 10, pos2 - pos1 - 10);
  double alpha = stod(sstr);

  ModelParameters parameters{alpha, frequencies, subst_params};
  return parameters;
}

} } // namespace pt::pll
