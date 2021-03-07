# Cuda NetCDF Interp

The aim of this tool is to accelerate 3D data interpolation using NVidia GPUs. For now I intend to implement a simple trilinear interpolate kernel, but more may follow.

## Installation

Install with CMake using install.sh

## Implementation Notes

In the present implementation, it'd be better to simply use a 1D block structure and compute points along the 1D array. However, for future kernels which might employ shared memory coalescence, it'd be better to have a 3D picture of the interpolation problem, so we accept some wastage at this stage.
## TODO

 - Implement trilinear interp kernel
 - Test templating for float/double
 - Check against e.g. numpy interp

