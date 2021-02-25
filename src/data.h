#include <string>
#include <netcdf.h>
#include <vector>
#include <map>
#include "types.h"

int getNetcdfData(std::string filename,
		    std::vector<std::string> dim_names,
		    std::vector<std::string> var_names,
		    std::vector<std::vector<double>> &coords,
		    std::vector<field_t<double>> &var_values);

gridspec_t getTestGrid();

gridspec_t getNetcdfGrid(coords_t ncCoords);

coords_t gridToGrid3D(gridspec_t inGrid, gridspec_t outGrid);

void cpuTrilinInterp(coords_t local_coords, std::vector<field_t<double>> fields, gridspec_t fieldsGrid);
