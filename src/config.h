#include <string>
#include "types.h"
#include <boost/program_options.hpp>

namespace po=boost::program_options;

config_t getCfg(std::string filename);

gridspec_t createOutGrid(const po::variables_map &varmap);
