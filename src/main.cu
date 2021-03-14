/* 3D NetCDF interpolation using CUDA

   Goal: Enable accelerated regridding of 3D data from NetCDF format

   Kernel: Trilinear interpolation w/ stored weights for repeat interp

   Challenges:

     - Off-axis grids
     - Calculation of block memory footprint
     - Negotiation of per-block parallel read

   Notes:

   For each point in dest_arr, only need to store id of single
   point in src_arr corresponding to bottom-left-near corner
   (min x,y,z) and <8> weights.

   If memory limited, could recompute 8 weights on-the-fly &
   only store <3> local coords.

   Each block corresponds to an (nx, ny, nz) subset of the dest_arr,
   and a corresponding subset of the src_arr (inc. ghost).

   Metis would take care of this, but let's do it manually as a
   learning exercise. Most efficient (minimize edge-cut)
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include <cuda_runtime_api.h>
#include <cuda.h>

#include "data.h"
#include "kernels.h"
#include "types.h"
#include "config.h"

using std::vector;
using std::cout;
using std::endl;

int main(int argc, char* argv[]){

  config_t config = getCfg(argv[1]);

  vector<std::string> dim_names = {"x","y","z"};
//   vector<std::string> var_names = {"temperature"};
//   // fields_t fields(var_names.size());
  vector<field_t<double>> fields(config.variables.size());
  coords_t coords(3);

  int retval = getNetcdfData(config.input_filename, config.dim_names, config.variables, coords, fields);
  if(retval){cout << "NetCDF Error" << endl; return retval;}

  // try the max function
  // ***********************************************
#if 0
  double *dfield, *dmax, *dtest;
  int nvals = fields[0].size();
  int blocksize = 512;
  int nblocks = (nvals + blocksize - 1) / blocksize;

  dtest = (double*)malloc(nvals*sizeof(double));
  gpuErrchk(cudaMalloc((void **)&dfield, nvals*sizeof(double)));
  gpuErrchk(cudaMalloc((void **)&dmax, nvals*sizeof(double))); // TODO too big!

  gpuErrchk(cudaMemcpy(dfield, &fields[0][0],nvals*sizeof(double), cudaMemcpyHostToDevice));

  // Repeated calls to kernel w/ block-level reduction
  int cnt = nvals;
  while(cnt > 1){
    cout << "Operating on cnt: " << cnt << endl;
    maxVal<<<nblocks, blocksize, blocksize*sizeof(double)>>>(cnt, dfield, dmax);
    gpuErrchk(cudaMemcpy(dfield, dmax, nvals*sizeof(double) ,cudaMemcpyDeviceToDevice));
    cnt = (cnt + blocksize - 1) / blocksize;
  }

  double maxVal;
  gpuErrchk(cudaMemcpy(&maxVal, dmax, sizeof(double), cudaMemcpyDeviceToHost));

  cout << "Max val from gpu: " << maxVal << endl;
  auto maxy = std::max_element(std::begin(fields[0]), std::end(fields[0]));
  cout << "Max val from cpu:" << *maxy << endl;

  // It works!
  // ***********************************************
#endif

  // Generate specification of new grid for interp (should be read in)
  gridspec_t gridSpecOut, gridSpecIn;
  gridSpecOut = config.gridspec_out;
  gridSpecIn = getNetcdfGrid(coords);

  // Copy fields to new vector & delete all values
  // Note that we thus ensure that dimension order is consistent
  // between input & output (fields[i].dim_order)
  vector<field_t<double>> interped_fields(fields.size()); // TODO proper copy constructor?
  for (unsigned int i = 0; i < fields.size(); i++){
    interped_fields[i].name = fields[i].name;
    interped_fields[i].dims = fields[i].dims;
    interped_fields[i].dim_order = fields[i].dim_order;
  }

  // CPU interp for benchmarking
  bool test = false;
  if(test){
    // Compute the coords of each point in new grid in the 'netcdf index' space
    coords_t localCoords = gridToGrid3D(gridSpecIn, gridSpecOut);
    cpuTrilinInterp(localCoords, fields, gridSpecIn, gridSpecOut, interped_fields);
  }

  gpuTrilinInterp(gridSpecIn, gridSpecOut, fields, interped_fields);

  retval = writeNetcdfData(config.output_filename, gridSpecOut, config.dim_names, interped_fields);
  if(retval){cout << "NetCDF Error" << endl; return retval;}

  return 0;

}
