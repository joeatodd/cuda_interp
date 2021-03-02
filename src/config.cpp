#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include "types.h"

namespace po=boost::program_options;


// Create the gridspec_t from input config
gridspec_t createOutGrid(const po::variables_map &varmap){
  gridspec_t outGridspec;

  outGridspec.nx[0] = varmap["Output.nx"].as<int>();
  outGridspec.nx[1] = varmap["Output.ny"].as<int>();
  outGridspec.nx[2] = varmap["Output.nz"].as<int>();

  outGridspec.dx[0] = varmap["Output.dx"].as<double>();
  outGridspec.dx[1] = varmap["Output.dy"].as<double>();
  outGridspec.dx[2] = varmap["Output.dz"].as<double>();

  outGridspec.x0[0] = varmap["Output.x0"].as<double>();
  outGridspec.x0[1] = varmap["Output.y0"].as<double>();
  outGridspec.x0[2] = varmap["Output.z0"].as<double>();

  return outGridspec;
}

// Get the configuration (filenames & desired grid)
config_t getCfg(std::string config_fname){
  config_t config;
  po::options_description opts("Options");

  opts.add_options() ("Input.filename", "Input filename");
  opts.add_options() ("Output.filename", "Output filename");
  opts.add_options() ("Output.x0", po::value<double>());
  opts.add_options() ("Output.y0", po::value<double>());
  opts.add_options() ("Output.z0", po::value<double>());

  opts.add_options() ("Output.dx", po::value<double>());
  opts.add_options() ("Output.dy", po::value<double>());
  opts.add_options() ("Output.dz", po::value<double>());

  opts.add_options() ("Output.nx", po::value<int>());
  opts.add_options() ("Output.ny", po::value<int>());
  opts.add_options() ("Output.nz", po::value<int>());


  std::ifstream config_stream(config_fname.c_str());
  po::variables_map vars;
  po::store(po::parse_config_file(config_stream, opts), vars);
  po::notify(vars);

  config.input_filename = vars["Input.filename"].as<std::string>();
  config.output_filename = vars["Output.filename"].as<std::string>();

  config.gridspec_out = createOutGrid(vars);

  std::cout << "Infile: " << vars["Input.filename"].as<std::string>() << std::endl;

  return config;
}
