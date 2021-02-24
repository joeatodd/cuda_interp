#include <string>
#include <netcdf.h>
#include <vector>

int get_netcdf_data(std::string filename,
		    std::vector<std::string> dim_names,
		    std::vector<std::string> var_names,
		    std::vector<std::vector<double>> &coords,
		    std::vector<std::vector<double>> &var_values);
