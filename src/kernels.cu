// CUDA kernels
#include "types.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <iostream>
#include <math.h>

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

// Max value (via recursive kernel calls)
// Makes use of recast trick because templated shared memory isn't allowed
template<typename T>
__global__
void maxVal(int N, T *arr, T *arr_out){
  extern __shared__ char smem[];
  T * shared = reinterpret_cast<T *>(smem);
  int gid = blockIdx.x*blockDim.x + threadIdx.x;
  int tid = threadIdx.x;

  // Load to shared mem then wait
  if(gid < N){
    shared[tid] = arr[gid];
  }
  __syncthreads();

  for(int i = (blockDim.x/2); i>0; i>>=1){
    if(tid < i && gid < N){
      T temp = shared[tid+i];
      if(temp > shared[tid]){
	shared[tid] = temp;
      }
    }
    __syncthreads();
  }

  if(gid < N && tid==0){
    arr_out[blockIdx.x] = shared[tid];
  }
}

__global__
void copyVal(int N, double *in_arr, double *out_arr){
  unsigned int gid = blockIdx.x*blockDim.x + threadIdx.x;
  if(gid < N){
    out_arr[gid] = in_arr[gid];
  }
}

/*
  First attempt at trilinear interpolation as CUDA kernel, using
  3D block structure
  inGrid - the grid specification of the input
  outGrid - the requested grid specification of the output
  di - order of dimensions (row-major)
*/
__global__
void trilin1(const gridspec_t *inGrid, const gridspec_t *outGrid, const int *di, double *inField, double *outField){

  // TODO - runs fail at the moment if all the code here is executed and any attempt is made to write
  // to outField[0] (or any other value?). So, memory corruption in the preceding code.
  // outField[0] = 1.0;
  unsigned int gid[3];
  gid[0] = blockIdx.x*blockDim.x + threadIdx.x;
  gid[1] = blockIdx.y*blockDim.y + threadIdx.y;
  gid[2] = blockIdx.z*blockDim.z + threadIdx.z;

  unsigned int nc[3];
  for (int i = 0; i < 3; i++) nc[i] = outGrid->nx[i];

  // Are we valid?
  bool active = gid[0] < nc[0] &&
    gid[1] < nc[1] &&
    gid[2] < nc[2];

  // Compute local coords & split into grid cells & in-cell coords
  double lc[3];
  for (int i = 0; i < 3; i++){
    double offset = outGrid->x0[i] - inGrid->x0[i];
    lc[i] = offset + (gid[di[i]] * outGrid->dx[di[i]])/inGrid->dx[di[i]];
  }

  double gridx_d[3];
  int gridx[3];
  for (int i=0; i<3; i++){
    lc[i] = modf(lc[i], &gridx_d[i]);
    gridx[i] = (int)gridx_d[i];
  }

  // Compute interpolation weights
  double weights[8];
  // C_000
  weights[0] = (1.0-lc[0])*(1.0-lc[1])*(1.0-lc[2]);

  // C_100
  weights[1] = (lc[0])*(1.0-lc[1])*(1.0-lc[2]);

  // C_010
  weights[2] = (1.0-lc[0])*(lc[1])*(1.0-lc[2]);

  // C_110
  weights[3] = (lc[0])*(lc[1])*(1-lc[2]);

  // C_001
  weights[4] = (1.0-lc[0])*(1.0-lc[1])*(lc[2]);

  // C_101
  weights[5] = (lc[0])*(1.0-lc[1])*(lc[2]);

  // C_011
  weights[6] = (1.0-lc[0])*(lc[1])*(lc[2]);

  // C_111
  weights[7] = lc[0]*lc[1]*lc[2];

  int idx[8];
  idx[0] = gridx[2]   + (gridx[1]   * nc[di[2]]) + ((gridx[0]   )* nc[di[1]] * nc[di[2]]);
  idx[1] = gridx[2]   + (gridx[1]   * nc[di[2]]) + ((gridx[0] +1)* nc[di[1]] * nc[di[2]]);
  idx[2] = gridx[2]   + (gridx[1]+1 * nc[di[2]]) + ((gridx[0]   )* nc[di[1]] * nc[di[2]]);
  idx[3] = gridx[2]   + (gridx[1]+1 * nc[di[2]]) + ((gridx[0] +1)* nc[di[1]] * nc[di[2]]);
  idx[4] = gridx[2]+1 + (gridx[1]   * nc[di[2]]) + ((gridx[0]   )* nc[di[1]] * nc[di[2]]);
  idx[5] = gridx[2]+1 + (gridx[1]   * nc[di[2]]) + ((gridx[0] +1)* nc[di[1]] * nc[di[2]]);
  idx[6] = gridx[2]+1 + (gridx[1]+1 * nc[di[2]]) + ((gridx[0]   )* nc[di[1]] * nc[di[2]]);
  idx[7] = gridx[2]+1 + (gridx[1]+1 * nc[di[2]]) + ((gridx[0] +1)* nc[di[1]] * nc[di[2]]);

  if(active){
    double value = 0.0;
    for (int n = 0; n < 8; n++){
      value += weights[n] * inField[idx[n]];
    }

    long unsigned int gcoord = gid[di[2]] + gid[di[1]]*nc[di[2]] + gid[di[0]] * nc[di[2]] * nc[di[1]];

    outField[gcoord] = value;
  }

}

__host__
void gpuTrilinInterp(const gridspec_t gridSpecIn, const gridspec_t gridSpecOut, const std::vector<field_t<double>> fieldsIn, std::vector<field_t<double>> fieldsOut){

  dim3 gridSize(2,2,2);
  dim3 blockSize(8,4,10);

  int gridMemsize = sizeof(gridspec_t);
  std::cout << "Gridspec_t size is: " << gridMemsize << std::endl;

  gridspec_t * dGridSpecIn, * dGridSpecOut;
  int * dDimOrder;

  gpuErrchk(cudaMalloc((void **)&dGridSpecIn, gridMemsize));
  gpuErrchk(cudaMalloc((void **)&dGridSpecOut, gridMemsize));
  gpuErrchk(cudaMalloc((void **)&dDimOrder, sizeof(int)*3));

  long unsigned int npts_out = gridSpecOut.nx[0] * gridSpecOut.nx[1] * gridSpecOut.nx[2];
  long unsigned int npts_in = gridSpecIn.nx[0] * gridSpecIn.nx[1] * gridSpecIn.nx[2];

  // // int npts = (2*2*2*8*4*10);
  // std::cout << "npts: " << npts << std::endl;

  double *dFieldOut, *dFieldIn;
  int inFieldMemSize = sizeof(double) * npts_in;
  int outFieldMemSize = sizeof(double) * npts_out;

  std::cout << "npts in: " << npts_in << " mem size: " << inFieldMemSize << std::endl;
  std::cout << "npts out: " << npts_out << " mem size: " << outFieldMemSize << std::endl;

  gpuErrchk(cudaMalloc((void **)&dFieldIn, inFieldMemSize));
  gpuErrchk(cudaMalloc((void **)&dFieldOut, outFieldMemSize));

  // localCoord = (double*)malloc(coordMemSize);

  // Copy data to GPU
  gpuErrchk(cudaMemcpy(dGridSpecIn, &gridSpecIn, gridMemsize, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dGridSpecOut, &gridSpecOut, gridMemsize, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dDimOrder, &fieldsIn[0].dim_order[0], sizeof(int)*3, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dFieldIn, &fieldsIn[0].values[0], inFieldMemSize, cudaMemcpyHostToDevice));

  trilin1<<<gridSize, blockSize>>>(dGridSpecIn, dGridSpecOut, &fieldsIn[0].dim_order[0], dFieldIn, dFieldOut);

  fieldsOut[0].values.resize(npts_out);
  gpuErrchk(cudaMemcpy(&fieldsOut[0].values[0], dFieldOut, outFieldMemSize, cudaMemcpyDeviceToHost));
}
