// module for loading 3D netcdf datasets

#include <iostream>
#include <string>
#include <netcdf.h>
#include <vector>
#include <assert.h>

#define netcdfERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

using std::vector;
using std::cout;
using std::endl;

/* Read a series of 3D arrays from a netcdf file
   This should be generalised to accept dimensions other than
   'x','y','z', and could also be used for data of other dimensionality.

   There is a NetCDF C++ interface but it doesn't seem to be available on
   the machine I'm using, so I'm mixing C++ and C code here.
*/
int get_netcdf_data(const std::string filename,
		    const vector<std::string> dim_names,
		    const vector<std::string> var_names,
		    vector<vector<double>> &coords,
		    vector<vector<double>> &var_values){

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

  // Read the dimension sizes (for now assume 'x','y','z')
  vector<size_t> dim_lens;
  vector<int> dim_ids;
  dim_lens.resize(3);
  dim_ids.resize(3);
  for (int i = 0; i < 3; i++){
    int dim_id;
    size_t dim_len; // TODO - directly pass ref to vector element here
    retval = nc_inq_dimid(ncid, dim_names[i].c_str(), &dim_id);
    if(retval) netcdfERR(retval);
    dim_ids[i] = dim_id;

    retval = nc_inq_dimlen(ncid, dim_ids[i], &dim_len);
    if(retval) netcdfERR(retval);
    dim_lens[i] = dim_len;
  }


  std::cout << "Input file has dimensions: ";
  int nvals = 1;
  for (int i = 0; i < 3; i++){
    std::cout << dim_lens[i] << " ";
    nvals *= dim_lens[i];
  }
  std::cout << endl;

  // TODO - generalise float/double
  // Read coordinates
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
  vector<int> var_ids(nvars_out);

  for (int i = 0; i < nvars_out; i++){
    int varid;
    retval = nc_inq_varid(ncid, var_names[i].c_str(), &varid);
    if(retval) netcdfERR(retval);

    var_values[i].resize(nvals);

    retval = nc_get_vara_double(ncid, varid, start, counts, &var_values[i][0]);
    if(retval) netcdfERR(retval);

  }

  return 0;
}
