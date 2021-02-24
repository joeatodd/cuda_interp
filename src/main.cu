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

#include "memops.cu"
#include "data.h"

using std::vector;
using std::cout;
using std::endl;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__global__
void multVal(int n, float coeff, float *x){
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n){
    x[idx] *= coeff; // In-place multiplication (is this allowed?)
  }
}

// 3D pitched equivalent of prev
__global__
void multVal3D(cudaPitchedPtr pitchPtr, cudaExtent ext, float coeff){
  size_t pitch = pitchPtr.pitch;


}

int main(){
  /* Declare dimensions for test array */
  int nx = 100, ny = 200, nz = 15; // Dimension
  double lx = 1000.0, ly = 1000.0, lz = 20.0; // Spacing
  double ox = 0.0, oy = 0.0, oz = 0.0; // Origin

  vector<std::string> dim_names = {"x","y","z"};
  vector<std::string> var_names = {"temperature"};
  vector<vector<double>> fields(var_names.size()), coords(3);

  int retval = get_netcdf_data("Seasonal_Spinup_8_body_temp.nc", dim_names, var_names, coords, fields);
  if(retval){cout << "NetCDF Error" << endl; return retval;}

  // TODO - determine dimension order here!

  /* Allocate the 3D test array */
  int N = nx * ny * nz;
  float *test_data;
  test_data = (float*)malloc(sizeof(float)*N); // TODO free

  // Try an actual 3D array rather than a flat one
  cudaExtent srcMemLayout = make_cudaExtent(sizeof(float)*nx, ny, nz);
  cudaPitchedPtr pitchPtr;
  gpuErrchk(cudaMalloc3D(&pitchPtr, srcMemLayout)); // TODO free

  std::cout << "Layout: " << srcMemLayout.width << " " << srcMemLayout.height << " " << srcMemLayout.depth << std::endl;

  // Set up copy params
  cudaMemcpy3DParms cp = {0};
  cp.srcPtr = make_cudaPitchedPtr(test_data, nx*sizeof(float), ny, nz);
  cp.dstPtr = pitchPtr;
  cp.extent = srcMemLayout;
  cp.srcPos = make_cudaPos(0,0,0);
  cp.dstPos = make_cudaPos(0,0,0);
  cp.kind = cudaMemcpyHostToDevice;

  // and equivalent for copy back
  cudaMemcpy3DParms pc = {0};
  cp.srcPtr = pitchPtr;
  cp.dstPtr = make_cudaPitchedPtr(test_data, nx*sizeof(float), ny, nz);
  cp.extent = srcMemLayout;
  cp.srcPos = make_cudaPos(0,0,0);
  cp.dstPos = make_cudaPos(0,0,0);
  cp.kind = cudaMemcpyDeviceToHost;


  /* Fill with 1s and copy to GPU */
  for (int i=0; i<N; i++){
    test_data[i] = 1.0;
  }

  // Copy!
  cudaMemcpy3D(&cp);

  // TODO compute occupancy
  float coeff = 2.5;
  /* multVal<<<(N+256-1)/256, 256>>>(N, coeff, pitchPtr); */
  multVal3D<<<(N+256-1)/256, 256>>>(pitchPtr, srcMemLayout, coeff);

  cudaMemcpy3D(&pc);
  /* cudaMemcpy(test_data, dtest_data, N*sizeof(float), cudaMemcpyDeviceToHost); */

  float err = 0.0;
  for (int i = 0; i < N; i++){

    err = max(err, abs(test_data[i] - coeff));
  }
  std::cout << "Max error: " << err << std::endl;
}
