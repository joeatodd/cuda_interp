# Cuda NetCDF Interp

The aim of this tool is to accelerate 3D data interpolation using NVidia GPUs. For now I intend to implement a simple trilinear interpolate kernel, but more may follow.

## Installation

Install with CMake using install.sh

## TODO

 - Define input format (or args?)
 - Implement trilinear interp kernel
 - Test templating for float/double
 - Check against e.g. numpy interp