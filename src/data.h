#include <string>
#include <netcdf.h>
#include <vector>
#include <map>
#include "types.h"

template <typename T>
int getNetcdfData(const std::string filename,
		  const std::vector<std::string> dim_names,
		  const std::vector<std::string> var_names,
		  std::vector<std::vector<T>> &coords,
		  std::vector<field_t<T>> &var_values);

template <typename T>
int writeNetcdfData(const std::string filename,
		     const gridspec_t gridspec,
		     const std::vector<std::string> dim_names,
		     const std::vector<field_t<T>> fields);

gridspec_t getTestGrid();

gridspec_t getNetcdfGrid(coords_t ncCoords);

coords_t gridToGrid3D(gridspec_t inGrid, gridspec_t outGrid);

template <typename T>
void cpuTrilinInterp(const coords_t local_coords,
		     const std::vector<field_t<T>> fields,
		     const gridspec_t fieldsGrid,
		     const gridspec_t outputGrid,
		     std::vector<field_t<T>> interped_fields);

int flatIdx(int i, int j, int k, int nj, int jk);
