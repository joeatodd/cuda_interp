// module for loading 3D netcdf datasets

#include <iostream>
#include <string>
#include <netcdf.h>
#include <vector>
#include <assert.h>
#include "types.h"

#define netcdfERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

using std::vector;
using std::cout;
using std::endl;

/* Read a series of 3D arrays from a netcdf file
   This should be generalised to accept dimensions other than
   'x','y','z', and could also be used for data of other dimensionality.

   There is a NetCDF C++ interface but it doesn't seem to be available on
   the machine I'm using, so I'm mixing C++ and C code here.

   @param filename NetCDF filename to read
   @param dim_names Prescribed dimensions of the data (TODO could read these instead)
   @param var_names Variable(s) to read
   @param coords The coordinates of the data (output)
   @param vars The data array for each variable
*/
int getNetcdfData(const std::string filename,
		  const vector<std::string> dim_names,
		  const vector<std::string> var_names,
		  vector<vector<double>> &coords,
		  vector<field_t<double>> &vars){

  int ncid, ndims, nvars_in, ngatts, unlimdim, nvars_out;
  int retval;

  nvars_out = var_names.size();

  // Open the NetCDF file & get basic info
  retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
  if(retval) netcdfERR(retval);
  retval = nc_inq(ncid, &ndims, &nvars_in, &ngatts, &unlimdim);
  if(retval) netcdfERR(retval);

  // Check we got 3D data & 3D dim list
  assert(dim_names.size() == 3);
  if (ndims != 3){
    printf("Error, expected 3 dimensions but found %i", ndims);
    return 2;
  }

  // Read the dimension sizes
  vector<size_t> dim_lens;
  vector<int> dim_ids;
  std::map<int, std::string> dim_map;

  dim_lens.resize(3);
  dim_ids.resize(3);
  for (int i = 0; i < 3; i++){
    int dim_id;
    size_t dim_len; // TODO - directly pass ref to vector element here
    retval = nc_inq_dimid(ncid, dim_names[i].c_str(), &dim_id);
    if(retval) netcdfERR(retval);
    dim_ids[i] = dim_id;
    dim_map[dim_id] = dim_names[i];

    retval = nc_inq_dimlen(ncid, dim_ids[i], &dim_len);
    if(retval) netcdfERR(retval);
    dim_lens[i] = dim_len;
  }

  // Report dimensions
  cout << "Input file has dimensions: " << endl;
  int nvals = 1;
  for (int i = 0; i < 3; i++){
    cout << dim_names[i] << " - ";
    cout << "id: " << dim_ids[i] << " ";
    cout << "size: " << dim_lens[i] << endl;
    nvals *= dim_lens[i];
  }
  cout << endl;

  // Read coordinates
  // TODO - generalise float/double
  for (int i = 0; i < 3; i++){
    coords[i].resize(dim_lens[i]);
    retval = nc_get_var_double(ncid, dim_ids[i], &coords[i][0]);
    if(retval) netcdfERR(retval);
  }

  // Start & counts (same for all vars)
  size_t start[3] = {0,0,0};
  size_t counts[3] = {dim_lens[2], dim_lens[1], dim_lens[0]};

  // Read the variables in turn (flattened)
  // vector<vector<double>> vars(nvars_out);
  for (int i = 0; i < nvars_out; i++){

    int varid;
    nc_type var_type;
    int varNDims;
    int varDims[NC_MAX_VAR_DIMS];

    // Get variable ID
    retval = nc_inq_varid(ncid, var_names[i].c_str(), &varid);
    if(retval) netcdfERR(retval);

    // Get the variable's dimensions
    nc_inq_var(ncid, varid, 0, &var_type, &varNDims, varDims, 0);
    assert(varNDims == 3);

    // Turn var_ids into names & store
    cout << "Variable " << var_names[i] << " has " << varNDims << " dimensions: ";
    vars[i].dims.resize(varNDims);
    for (int j = 0; j < varNDims; j++){
      vars[i].dims[j] = dim_map[varDims[j]];
      cout << vars[i].dims[j] << " ";
    }
    cout << endl;

    vars[i].values.resize(nvals);

    // Read in the values
    retval = nc_get_vara_double(ncid, varid, start, counts, &vars[i].values[0]);
    if(retval) netcdfERR(retval);
  }

  return 0;
}

// Create an example grid spec for 3D interpolation
gridspec_t getTestGrid(){

  gridspec_t my_gridspec;

  my_gridspec.dx[0] = 50.0;
  my_gridspec.dx[1] = 50.0;
  my_gridspec.dx[2] = 10.0;

  my_gridspec.x0[0] = -210000.0;
  my_gridspec.x0[1] = -2135000.0;
  my_gridspec.x0[2] = -1000.0;

  my_gridspec.nx[0] = 1000.0;
  my_gridspec.nx[1] = 1000.0;
  my_gridspec.nx[2] = 100.0;

  return my_gridspec;
}

// Given a vector<vector<...>> of coordinates, produce the gridspec
gridspec_t getNetcdfGrid(coords_t coords){

  gridspec_t ncGridspec;

  // X, Y, Z
  for (int i = 0; i < 3; i++){
    ncGridspec.x0[i] = coords[i][0];
    ncGridspec.dx[i] = coords[i][1] - coords[i][0];
    ncGridspec.nx[i] = coords[i].size();
  }

  return ncGridspec;
}

/*
Compute the coordinates of each output point in the
netCDF index space. So, a point with local_coords = 2.5, 30.7, 19.7
sits between the 8 input points with x = 2,3,  y=30,31, z=19,20

Note - generalisation to rotated coordinates would require full
specification of all 3 coordinates for every single output point.
*/
coords_t gridToGrid3D(gridspec_t inGrid, gridspec_t outGrid){

  cout << "Hello from gridtoGrid3D" << endl;
  coords_t local_coords(3);

  double offsets[3];

  for (int i = 0; i < 3; i++){
    // Compute origin
    offsets[i] = outGrid.x0[i] - inGrid.x0[i];
    // pre-allocate space
    local_coords[i].resize(outGrid.nx[i]);
    for (int j = 0; j < outGrid.nx[i]; j++){
      local_coords[i][j] = (offsets[i] + j * outGrid.dx[i]) / inGrid.dx[i];
    }
  }

  return local_coords;
}

/*
  Compute trilinear interpolation of fields in CPU
  @param local_coords: coordinates of output points in 'input grid' index coord system
  @param fields: the input data
  @param fieldsGrid: the gridspec_t object defining origin, nx, dx of the input

*/
void cpuTrilinInterp(coords_t local_coords, vector<field_t<double>> fields, gridspec_t fieldsGrid){
  int nFields = fields.size();

}
