#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <map>
#include <vector>

struct gridspec_t{
  int nx[3];
  double dx[3], x0[3];
};

template <typename T>
struct field_t{
  std::vector<T> values;
  std::vector<std::string> dims;
};

typedef std::vector<std::vector<double>> coords_t;
typedef std::vector<std::vector<double>> fields_t;

#endif
