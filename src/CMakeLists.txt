add_executable(CudaInterp main.cu data.cpp kernels.cu config.cpp)

target_include_directories(CudaInterp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(NetCDF REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

include_directories( ${Boost_INCLUDE_DIRS} )
target_link_libraries(CudaInterp netcdf ${Boost_LIBRARIES})

target_compile_options(CudaInterp PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_70,code=sm_70>)
