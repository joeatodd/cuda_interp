__global__
void multVal(int n, float coeff, float *x);

__global__
void multVal3D(cudaPitchedPtr pitchPtr, cudaExtent ext, float coeff);

template<typename T>
__global__
void maxVal(int N, T *arr, T *arr_out);

__global__
void copyVal(int N, double *in_arr, double *out_arr);
