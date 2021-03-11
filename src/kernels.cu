// CUDA kernels
#include "types.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <assert.h>

using std::vector;
using std::cout;
using std::endl;

__global__
void multVal(int n, float coeff, float *x){
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n){
    x[idx] *= coeff; // In-place multiplication (is this allowed?)
  }
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

/*
  First attempt at trilinear interpolation as CUDA kernel, using
  3D block structure
  inGrid - the grid specification of the input
  outGrid - the requested grid specification of the output
  di - order of dimensions (row-major)
*/
__global__
void trilin1(const gridspec_t *inGrid, const gridspec_t *outGrid, const int *di, double *inField, double *outField){

  // 3D global index
  unsigned int gid[3];
  gid[0] = blockIdx.x*blockDim.x + threadIdx.x;
  gid[1] = blockIdx.y*blockDim.y + threadIdx.y;
  gid[2] = blockIdx.z*blockDim.z + threadIdx.z;

  // Copy grid specs to registers
  unsigned int nc_in[3], nc_out[3];
  for (int i = 0; i < 3; i++){
    nc_out[i] = outGrid->nx[i];
    nc_in[i]  = inGrid->nx[i];
  }

  // Are we valid?
  bool active = gid[0] < nc_out[0] &&
    gid[1] < nc_out[1] &&
    gid[2] < nc_out[2];

  // Compute local coords & split into grid cells & in-cell coords
  // Here lc is ordered X, Y, Z
  double lc[3];
  for (int i = 0; i < 3; i++){
    double offset = outGrid->x0[i] - inGrid->x0[i];
    lc[i] = (offset + gid[i] * outGrid->dx[i])/inGrid->dx[i];
  }

  // Compute surrounding inGrid indices & grid-local coords
  // Here ordered X, Y, Z (inherit from lc)
  double gridx_d[3];
  int gridx[3];
  for (int i=0; i<3; i++){
    lc[i] = modf(lc[i], &gridx_d[i]);
    gridx[i] = static_cast<int>(gridx_d[i]);
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


  // Weights are computed in X, Y, Z (not Z, Y, X), so we need to get same order idx
  // idx[0] is x0, y0, z0
  // idx[1] is x1, y0, z0
  int idx[8];
  idx[0] = gridx[di[2]]   + ((gridx[di[1]]  ) * nc_in[di[2]]) + ((gridx[di[0]]   )* nc_in[di[1]] * nc_in[di[2]]);
  idx[1] = gridx[di[2]]+1 + ((gridx[di[1]]  ) * nc_in[di[2]]) + ((gridx[di[0]]   )* nc_in[di[1]] * nc_in[di[2]]);
  idx[2] = gridx[di[2]]   + ((gridx[di[1]]+1) * nc_in[di[2]]) + ((gridx[di[0]]   )* nc_in[di[1]] * nc_in[di[2]]);
  idx[3] = gridx[di[2]]+1 + ((gridx[di[1]]+1) * nc_in[di[2]]) + ((gridx[di[0]]   )* nc_in[di[1]] * nc_in[di[2]]);
  idx[4] = gridx[di[2]]   + ((gridx[di[1]]  ) * nc_in[di[2]]) + ((gridx[di[0]] +1)* nc_in[di[1]] * nc_in[di[2]]);
  idx[5] = gridx[di[2]]+1 + ((gridx[di[1]]  ) * nc_in[di[2]]) + ((gridx[di[0]] +1)* nc_in[di[1]] * nc_in[di[2]]);
  idx[6] = gridx[di[2]]   + ((gridx[di[1]]+1) * nc_in[di[2]]) + ((gridx[di[0]] +1)* nc_in[di[1]] * nc_in[di[2]]);
  idx[7] = gridx[di[2]]+1 + ((gridx[di[1]]+1) * nc_in[di[2]]) + ((gridx[di[0]] +1)* nc_in[di[1]] * nc_in[di[2]]);

  double weightSum = 0.0;
  if(active){
    double value = 0.0;
    for (int n = 0; n < 8; n++){
      value += weights[n] * inField[idx[n]];
      weightSum += weights[n];
    }
    assert(fabs(weightSum - 1) < 1.0e-10);
    // Compute flat idx for output
    unsigned int gcoord = gid[di[2]] + gid[di[1]]*nc_out[di[2]] +
      gid[di[0]] * nc_out[di[2]] * nc_out[di[1]];
    outField[gcoord] = value;
  }

}

/*
Computes the (naturally) ordered list of prime factors of n
*/
vector<int> primeFactors(int n){
  vector<int> factors;

  int root_n = sqrt(n);

  while(n % 2 == 0){
    n /= 2;
    factors.push_back(2);
  }

  for (int i = 3; i <= root_n; i+=2){
    while(n % i == 0){
      n /= i;
      factors.push_back(i);
    }
  }

  if(n > 1) factors.push_back(n);
  return factors;
}


/*
Compute the optimal block & grid dimensions to 1) maximise occupancy on device and
2) minimise inactive threads. This is a tricky optimisation problem involving prime
factorisation.

First compute the prime factors of the blockSize (lots of 2s!), then find which
dimension of the gridspec_t best fits this factor (least waste in other 2 dims)
*/
__host__
void computeBlockGridDims(const gridspec_t &gridSpecOut, const int blockSizeOcc, dim3 &gridDims, dim3 &blockDims){

  vector<int> blockFactors = primeFactors(blockSizeOcc);

  std::array<int,3> currDims = {1,1,1};
  std::array<int,3> gridspecDims = {gridSpecOut.nx[0], gridSpecOut.nx[1], gridSpecOut.nx[2]};

  // Cycle prime factors of block size
  for (auto factor : blockFactors){

    std::array<int,3> remainders;
    for (int i = 0; i < 3; i++){
      int n = currDims[i] * factor;
      remainders[i] = n - (gridspecDims[i] % n); // how many rows left empty at edge?
    }

    // For each possible remainder, the total wasted threads is the
    // product of the remainder and the other two dimensions of the grid
    std::array<int,3> loss;
    for (int i = 0; i < 3; i++){
      loss[i] = remainders[i];
      for (int j = 0; j < 3; j++){
	if(i==j) continue;
	loss[i] *= gridspecDims[j];
      }
    }
    // Which fits best?
    int bestIdx = std::distance(loss.begin(), std::min_element(loss.begin(), loss.end()));

    currDims[bestIdx] *= factor;
  }

  std::array<double, 3> gridDimsFloat;
  std::array<int,3> gridDimsArr;
  int totalBlocks = 1;
  for (int i = 0; i < 3; i++){
    gridDimsFloat[i] = static_cast<double>(gridspecDims[i]) / static_cast<double>(currDims[i]);
    gridDimsArr[i] = static_cast<int>(ceil(gridDimsFloat[i]));
    totalBlocks *= gridDimsArr[i];
  }

  blockDims.x = currDims[0];
  blockDims.y = currDims[1];
  blockDims.z = currDims[2];

  gridDims.x = gridDimsArr[0];
  gridDims.y = gridDimsArr[1];
  gridDims.z = gridDimsArr[2];

  int nPts = gridspecDims[0] * gridspecDims[1] * gridspecDims[2];
  int theoryBlocks = (nPts + blockSizeOcc - 1) / blockSizeOcc;

  cout << "Theoretical min blocks: " << theoryBlocks << endl;
  cout << "Actual num blocks: "  << totalBlocks << endl;
  cout << "Grid efficiency: " << static_cast<double>(theoryBlocks) / static_cast<double>(totalBlocks);
  cout << endl;
}

__host__
void gpuTrilinInterp(const gridspec_t &gridSpecIn, const gridspec_t &gridSpecOut, const std::vector<field_t<double>> &fieldsIn, std::vector<field_t<double>> &fieldsOut){

  int grid_tMemSize = sizeof(gridspec_t);
  std::cout << "Gridspec_t size is: " << grid_tMemSize << std::endl;

  gridspec_t * dGridSpecIn, * dGridSpecOut;
  int * dDimOrder;

  gpuErrchk(cudaMalloc((void **)&dGridSpecIn, grid_tMemSize));
  gpuErrchk(cudaMalloc((void **)&dGridSpecOut, grid_tMemSize));
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

  std::cout << "GridSpecIn nx: " << gridSpecIn.nx[0] << " " << gridSpecIn.nx[1] << " " << gridSpecIn.nx[2] << std::endl;
  std::cout << "GridSpecOut nx: " << gridSpecOut.nx[0] << " " << gridSpecOut.nx[1] << " " << gridSpecOut.nx[2] << std::endl;

  std::cout << "GridSpecIn dx: " << gridSpecIn.dx[0] << " " << gridSpecIn.dx[1] << " " << gridSpecIn.dx[2] << std::endl;
  std::cout << "GridSpecOut dx: " << gridSpecOut.dx[0] << " " << gridSpecOut.dx[1] << " " << gridSpecOut.dx[2] << std::endl;

  gpuErrchk(cudaMalloc((void **)&dFieldIn, inFieldMemSize));
  gpuErrchk(cudaMalloc((void **)&dFieldOut, outFieldMemSize));

  // localCoord = (double*)malloc(coordMemSize);

  // Copy data to GPU
  gpuErrchk(cudaMemcpy(dGridSpecIn, &gridSpecIn, grid_tMemSize, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dGridSpecOut, &gridSpecOut, grid_tMemSize, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dDimOrder, &fieldsIn[0].dim_order[0], sizeof(int)*3, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(dFieldIn, &fieldsIn[0].values[0], inFieldMemSize, cudaMemcpyHostToDevice));

  // Block size calculator
  int blockSize1D, minGridSize;
  cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize1D, trilin1, 0, npts_out);

  std::cout << "Occupancy results, blocksize: " << blockSize1D << std::endl;

  dim3 gridSize, blockSize;
  computeBlockGridDims(gridSpecOut, blockSize1D, gridSize, blockSize);

  std::cout << "Grid size: " << gridSize.x << " " <<  gridSize.y << " " <<  gridSize.z << std::endl;
  std::cout << "Block size: " << blockSize.x << " " <<  blockSize.y << " " <<  blockSize.z << std::endl;

  trilin1<<<gridSize, blockSize>>>(dGridSpecIn, dGridSpecOut, dDimOrder, dFieldIn, dFieldOut);

  fieldsOut[0].values.resize(npts_out);

  gpuErrchk(cudaMemcpy(&fieldsOut[0].values[0], dFieldOut, outFieldMemSize, cudaMemcpyDeviceToHost));
}
