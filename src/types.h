#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <map>
#include <vector>

struct gridspec_t{
  std::vector<int> nx;
  std::vector<double> dx, x0;
};

template <typename T>
struct field_t{
  std::string name;
  std::vector<T> values;
  std::vector<std::string> dims;
  std::vector<int> dim_order;
};

typedef std::vector<std::vector<double>> coords_t;

#endif
