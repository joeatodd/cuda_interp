nvcc -gencode arch=compute_70,code=sm_70 ./src/main.cu ./src/data.cpp -lnetcdf -o interp.so
