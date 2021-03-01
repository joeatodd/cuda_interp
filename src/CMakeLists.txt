add_executable(CudaInterp main.cu data.cpp kernels.cu)

target_include_directories(CudaInterp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(CudaInterp netcdf)

target_compile_options(CudaInterp PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_70,code=sm_70>)
