# Cuda NetCDF Interp

The aim of this tool is to accelerate 3D data interpolation using NVidia GPUs. For now, the kernel is just a simple trilinear interpolation, but more may follow.

## Installation

Build & install with CMake using install.sh. By default this installs the project into build/bin; this is controlled via the DCMAKE\_INSTALL\_PREFIX argument.

## Implementation Notes

### Grid & Block Structure

As CudaInterp addresses 3D interpolation problems, it adopts a 3D CUDA block/grid structure. In the present implementation, it'd actually be more efficient to simply use a 1D block structure and compute points along the 1D array. However, for future kernels which will employ shared memory coalescence for loading input data, it's be better to have a 3D picture of the interpolation problem, so we accept some wastage at this stage.

We make use of the Cuda Occupancy API to query the best (scalar) blocksize for the target device, then estimate the best 3D block dimensions via prime factorisation.

### NetCDF Interface

The NetCDF C interface (as opposed to C++) is used, as this seems to be more widely available.

## TODO

 - Fix templating for float/double
 - Check against e.g. numpy interp
 - Write interfaces
 - Other interpolation algorithms (e.g. tricubic)
 - Diagonal grids & polar coordinates
