// module for loading 3D netcdf datasets

#include <iostream>
#include <string>
#include <netcdf.h>
#include <vector>
#include <assert.h>
#include "types.h"
#include <math.h>
#include <numeric>

#define netcdfERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

using std::vector;
using std::cout;
using std::endl;

// Returns a 3D -> flat array index assuming C ordering (last varies fastest)
int flatIdx(int i, int j, int k, int nj, int nk){
  // result = i + (j * ni) + (k * ni * nj);
  return k + (j * nk) + (i * nj * nk);
}

/* Read a series of 3D arrays from a netcdf file
   This should be generalised to accept dimensions other than
   'x','y','z', and could also be used for data of other dimensionality.

   There is a NetCDF C++ interface but it doesn't seem to be available on
   the machine I'm using, so I'm mixing C++ and C code here.

   input  filename NetCDF filename to read
   input  dim_names Prescribed dimensions of the data (TODO could read these instead)
   input  var_names Variable(s) to read
   output  coords The coordinates of the data (output)
   output  vars The data array for each variable
*/
template <typename T>
int getNetcdfData(const std::string filename,
		  const vector<std::string> dim_names,
		  const vector<std::string> var_names,
		  vector<vector<T>> &coords,
		  vector<field_t<T>> &vars){

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
    retval = nc_inq_dimid(ncid, dim_names[i].c_str(), &dim_ids[i]);
    if(retval) netcdfERR(retval);
    dim_map[dim_ids[i]] = dim_names[i];

    retval = nc_inq_dimlen(ncid, dim_ids[i], &dim_lens[i]);
    if(retval) netcdfERR(retval);
  }

  // Report (& check) dimensions
  cout << "Input file has dimensions: " << endl;
  int nvals = 1;
  for (int i = 0; i < 3; i++){
    cout << dim_names[i] << " - ";
    cout << "id: " << dim_ids[i] << " ";
    cout << "size: " << dim_lens[i] << endl;

    nvals *= dim_lens[i];
    assert(dim_ids[i] == i); // NetCDF docs specify dims are numbered from 0 in C.
  }
  cout << endl;

  // Read coordinates
  for (int i = 0; i < 3; i++){
    coords[i].resize(dim_lens[i]);
    retval = nc_get_var(ncid, dim_ids[i], &coords[i][0]);
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

    vars[i].name = var_names[i];

    // Get variable ID
    retval = nc_inq_varid(ncid, var_names[i].c_str(), &varid);
    if(retval) netcdfERR(retval);

    // Get the variable's dimensions
    nc_inq_var(ncid, varid, 0, &var_type, &varNDims, varDims, 0);
    assert(varNDims == 3);

    // Turn var_ids into names & store
    vars[i].dims.resize(varNDims);
    vars[i].dim_order.resize(varNDims);
    for (int j = 0; j < varNDims; j++){
      vars[i].dim_order[j] = varDims[j];
      vars[i].dims[j] = dim_map[varDims[j]];
    }

    // Print out variable dim info
    cout << "Variable " << var_names[i] << " has " << varNDims << " dimensions: ";
    for (int j = 0; j < varNDims; j++){
      cout << vars[i].dims[j] << "(" << vars[i].dim_order[j] << ") ";
    }
    cout << endl;

    // Read in the values
    vars[i].values.resize(nvals);
    retval = nc_get_vara(ncid, varid, start, counts, &vars[i].values[0]);
    if(retval) netcdfERR(retval);
  }

  return 0;
}


// TODO - this is a template function but only writes doubles...
template <typename T>
int writeNetcdfData(const std::string filename,
		     const gridspec_t gridspec,
		     const vector<std::string> dim_names,
		     const vector<field_t<T>> fields){

  int nDims = dim_names.size();
  int nVars = fields.size();

  // Create the file
  int ncid, retval;
  retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  if(retval) netcdfERR(retval);

  // Define dims & coord vars
  vector<int> dim_ids(nDims), dim_var_ids(nDims);
  for (int i = 0; i < nDims; i++){
    // Define as dimension
    retval = nc_def_dim(ncid, dim_names[i].c_str(), gridspec.nx[i], &dim_ids[i]);
    if(retval) netcdfERR(retval);

    // Define as variable
    retval = nc_def_var(ncid, dim_names[i].c_str(), NC_DOUBLE, 1, &dim_ids[i], &dim_var_ids[i]);
    if(retval) netcdfERR(retval);
  }
  // TODO - write attributes to dimensions?

  // Define field vars
  vector<int> var_ids(nVars);
  for (int i = 0; i < nVars; i++){

    // Create dim permutation for this variable
    vector<int> var_dim_ids(nDims);
    for (int j = 0; j < nDims; j++) var_dim_ids[j] = dim_ids[fields[i].dim_order[j]];
    retval = nc_def_var(ncid, fields[i].name.c_str(), NC_DOUBLE, nDims, &var_dim_ids[0], &var_ids[i]);
    if(retval) netcdfERR(retval);
  }

  // End 'define' mode
  retval = nc_enddef(ncid);
  if(retval) netcdfERR(retval);

  // Write dims
  for (int i = 0; i < nDims; i++){
    vector<double> coord_vals(gridspec.nx[i]);

    for (int j = 0; j < gridspec.nx[i]; j++)
      coord_vals[j] = gridspec.x0[i] + j*gridspec.dx[i];

    retval = nc_put_var_double(ncid, dim_var_ids[i], &coord_vals[0]);
    if(retval) netcdfERR(retval);
  }

  // Write vars
  for (int i = 0; i < nVars; i++){
    cout << "Debug first value: " << fields[i].values[0] << endl;
    retval = nc_put_var_double(ncid, var_ids[i], &fields[i].values[0]);
    if(retval) netcdfERR(retval);
  }

  retval = nc_close(ncid);
  if(retval) netcdfERR(retval);

  return 0;
}

// Create an example grid spec for 3D interpolation
gridspec_t getTestGrid(){

  gridspec_t my_gridspec;

  my_gridspec.nx.resize(3);
  my_gridspec.dx.resize(3);
  my_gridspec.x0.resize(3);

  my_gridspec.dx[0] = 50.0;
  my_gridspec.dx[1] = 50.0;
  my_gridspec.dx[2] = 10.0;

  my_gridspec.x0[0] = -210000.0;
  my_gridspec.x0[1] = -2135000.0;
  my_gridspec.x0[2] = -200.0;

  my_gridspec.nx[0] = 150;
  my_gridspec.nx[1] = 100;
  my_gridspec.nx[2] = 10;

  // 1000, 1000, 100 (i.e. x, y, z)
  cout << "Test Grid dims: ";
  for (int i = 0; i < 3; i++){
    cout << my_gridspec.nx[i] << " ";
  }
  cout << endl;

  return my_gridspec;
}

// Given a vector<vector<...>> of coordinates, produce the gridspec
gridspec_t getNetcdfGrid(coords_t coords){

  gridspec_t ncGridspec;

  ncGridspec.nx.resize(3);
  ncGridspec.dx.resize(3);
  ncGridspec.x0.resize(3);

  // X, Y, Z
  for (int i = 0; i < 3; i++){
    ncGridspec.x0[i] = coords[i][0];
    ncGridspec.dx[i] = coords[i][1] - coords[i][0];
    ncGridspec.nx[i] = coords[i].size();
  }

  // X, Y, Z (i.e. not as ordered in variable, as defined in file)
  cout << "NetCDF Grid dims: ";
  for (int i = 0; i < 3; i++){
    cout << ncGridspec.nx[i] << " ";
  }
  cout << endl;

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

  Each new point is interpolated from 8 surrounding points which we can label C_xyz:

     C_000, C_100, C_010, C_110, C_001, C_101, C_011, C_111

  @param local_coords: coordinates of output points in 'input grid' index coord system
  @param fields: the input data
  @param fieldsGrid: the gridspec_t object defining origin, nx, dx of the input

*/
template <typename T>
void cpuTrilinInterp(const coords_t local_coords, const vector<field_t<T>> fields,
		     const gridspec_t fieldsGrid, const gridspec_t outputGrid,
		     vector<field_t<T>> interped_fields){

  // How many fields, dimensions and values?
  int nFields = fields.size();
  int nDims = local_coords.size();
  int nVals = 1;
  for (int i = 0; i < nDims; i++){
    nVals *= local_coords[i].size();
  }

  // Alias di for dimension order
  const vector<int>& di = fields[0].dim_order;

  // Preallocate space for values
  for (int i = 0; i < nFields; i++){
    interped_fields[i].values.resize(nVals);
  }

  cout << "Local coords dims: ";
  for (int i = 0; i < 3; i++){
    cout << i << ": " << local_coords[i].size() << ", ";
  }
  cout << endl;

  cout << "di: ";
  for (int i = 0; i < 3; i++){
    cout << di[i] << " ";
  }
  cout << endl;

  cout << "fieldsGrid nx: ";
  for (int i = 0; i < 3; i++){
    cout << fieldsGrid.nx[i] << " ";
  }
  cout << endl;

  cout << "outputGrid nx: ";
  for (int i = 0; i < 3; i++){
    cout << outputGrid.nx[i] << " ";
  }
  cout << endl;

  cout << "Local coords range: ";
  for (int i = 0; i < 3; i++){
    cout << local_coords[i][0] << " ";
    cout << local_coords[i][local_coords[i].size()-1] << ", ";
  }
  cout << endl;

  T x, y, z, weights[8];
  T dii, djj, dkk;
  int ii, jj, kk, idx[8];

  bool printed = false;

  for(unsigned int i = 0; i < local_coords[di[0]].size(); i++){
    for(unsigned int j = 0; j < local_coords[di[1]].size(); j++){
      for(unsigned int k = 0; k < local_coords[di[2]].size(); k++){


	// Split local coord into integer index & in-box float
	x = modf(local_coords[di[0]][i], &dii);
	y = modf(local_coords[di[1]][j], &djj);
	z = modf(local_coords[di[2]][k], &dkk);
	ii = (int)dii; jj = (int)djj; kk = (int)dkk;

	// Check local coord within bounds
	if (ii < 0 || ii >= fieldsGrid.nx[di[0]] - 1 ||
	    jj < 0 || jj >= fieldsGrid.nx[di[1]] - 1 ||
	    kk < 0 || kk >= fieldsGrid.nx[di[2]] - 1) continue;

	// C_000
	weights[0] = (1.0-x)*(1.0-y)*(1.0-z);

	// C_100
	weights[1] = (x)*(1.0-y)*(1.0-z);

	// C_010
	weights[2] = (1.0-x)*(y)*(1.0-z);

	// C_110
	weights[3] = (x)*(y)*(1-z);

	// C_001
	weights[4] = (1.0-x)*(1.0-y)*(z);

	// C_101
	weights[5] = (x)*(1.0-y)*(z);

	// C_011
	weights[6] = (1.0-x)*(y)*(z);

	// C_111
	weights[7] = x*y*z;

	T weightsum = 0;
	weightsum = std::accumulate(weights, weights+8, weightsum);
	assert(fabs(weightsum - 1) < 1.0e-10); // TODO get actual machine eps

	idx[0] = flatIdx(ii,   jj,   kk,   fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[1] = flatIdx(ii+1, jj,   kk,   fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[2] = flatIdx(ii,   jj+1, kk,   fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[3] = flatIdx(ii+1, jj+1, kk,   fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[4] = flatIdx(ii,   jj,   kk+1, fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[5] = flatIdx(ii+1, jj,   kk+1, fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[6] = flatIdx(ii,   jj+1, kk+1, fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);
	idx[7] = flatIdx(ii+1, jj+1, kk+1, fieldsGrid.nx[di[1]], fieldsGrid.nx[di[2]]);

	if(!printed){
	  printed = true;
	  cout << "idx: ";
	  for (int z = 0; z < 8; z++) cout << idx[z] <<  " ";
	  cout << endl;
	}

	// Which index of the output field are we filling?
	int outidx = flatIdx(i,j,k, outputGrid.nx[di[1]], outputGrid.nx[di[2]]);

	// Interp each field in turn (note, this relies on all fields defined on same dims)
	for (int f = 0; f < nFields; f++){
	  T value = 0.0;
	  for (int n = 0; n < 8; n++){
	    value += weights[n] * fields[f].values[idx[n]];
	  }
	  interped_fields[f].values[outidx] = value;
	  // cout << "Writing value " << value << " to idx " << outidx << endl;
	}
      }
    }
  }
}


// Template instantiation
template int getNetcdfData(const std::string filename,
			   const std::vector<std::string> dim_names,
			   const std::vector<std::string> var_names,
			   std::vector<std::vector<double>> &coords,
			   std::vector<field_t<double>> &var_values);

template int getNetcdfData(const std::string filename,
			   const std::vector<std::string> dim_names,
			   const std::vector<std::string> var_names,
			   std::vector<std::vector<float>> &coords,
			   std::vector<field_t<float>> &var_values);

template void cpuTrilinInterp(const coords_t local_coords,
			      const std::vector<field_t<float>> fields,
			      const gridspec_t fieldsGrid,
			      const gridspec_t outGrid,
			      const std::vector<field_t<float>> interped_fields);

template void cpuTrilinInterp(const coords_t local_coords,
			      const std::vector<field_t<double>> fields,
			      const gridspec_t fieldsGrid,
			      const gridspec_t outGrid,
			      std::vector<field_t<double>> interped_fields);

template int writeNetcdfData(const std::string filename,
		     const gridspec_t gridspec,
		     const vector<std::string> dim_names,
		     const vector<field_t<double>> fields);
