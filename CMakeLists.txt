cmake_minimum_required(VERSION 3.9)

# set the project name
project(CudaInterp LANGUAGES CUDA C CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Location of FindNetCDF
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMakeFiles")

find_package(NetCDF REQUIRED)

# Include the source files
add_subdirectory (src)
