// CUDA kernels

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
