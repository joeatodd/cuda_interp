#include "types.h"

__global__
void multVal(int n, float coeff, float *x);

__global__
void multVal3D(cudaPitchedPtr pitchPtr, cudaExtent ext, float coeff);

template<typename T>
__global__
void maxVal(int N, T *arr, T *arr_out);

__global__
void copyVal(int N, double *in_arr, double *out_arr);

__global__
void trilin1(gridspec_t *inGrid, gridspec_t *outGrid);

__host__
void gpuTrilinInterp(const gridspec_t &gridSpecIn, const gridspec_t &gridspecOut, const std::vector<field_t<double>> &fieldsIn, std::vector<field_t<double>> &fieldsOut);
