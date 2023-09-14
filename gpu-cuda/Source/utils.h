// Host code
__host__ void PrintGPUmem();
__host__ int getSPcores(cudaDeviceProp devProp);
__host__ void throwWARNING(const std::string messg);
__host__ void throwERROR(const std::string messg);
