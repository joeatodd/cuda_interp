#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <map>
#include <vector>
#include <cuda_runtime_api.h>
#include <cuda.h>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


struct gridspec_t{
  int nx[3];
  double dx[3], x0[3];
};

template <typename T>
struct field_t{
  std::string name;
  std::vector<T> values;
  std::vector<std::string> dims;
  std::vector<int> dim_order;
};

struct config_t{
  std::string input_filename;
  std::string output_filename;
  gridspec_t gridspec_out;
};

typedef std::vector<std::vector<double>> coords_t;

#endif
