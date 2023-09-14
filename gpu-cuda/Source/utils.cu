//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// Utility functions
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains utility functions for CPU/GPU code
// PST, last update 05/30/2017

#include <iostream>
#include <string>
#include <exception>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "utils.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-specific functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void PrintGPUmem()

////////////////////////////////////////////////////////////////////////
// Prints device memory stats

{
// Memory info
   size_t freem, totalm;
   prec_typ free_m,total_m,used_m;
//   cudaMemGetInfo((size_t*)&freem,(size_t*)&totalm);
   cudaMemGetInfo(&freem,&totalm);
   free_m =(size_t)freem/1048576.0;
   total_m=(size_t)totalm/1048576.0;
   used_m=(total_m-free_m);
   printf(" device total memory: %15.2f MB \n",total_m);
   printf(" device  used memory: %15.2f MB \n",used_m);
   printf(" device  free memory: %15.2f MB \n",free_m);

// Cores info
   int regs,cores,nDevices;
   cudaDeviceProp devProp;
   cudaGetDeviceCount(&nDevices);
   printf(" *** %d devices detected *** \n",nDevices);
   for (int i=0;i<nDevices;i++){
       cudaGetDeviceProperties(&devProp,i);
       cores=getSPcores(devProp);
       regs=devProp.regsPerBlock/BLOCK_SIZE_2D/BLOCK_SIZE_2D;
       printf(" device *** %d *** shared memory per block: %lu\n",i,devProp.sharedMemPerBlock);
       printf(" device *** %d *** registers per thread:      %d \n",i,regs);
       printf(" device *** %d *** streaming multiprocessors: %d \n",i,devProp.multiProcessorCount);
       printf(" device *** %d *** cores: %d \n",i,cores);
   }
}

////////////////////////////////////////////////////////////////////////

__host__ int getSPcores(cudaDeviceProp devProp)

////////////////////////////////////////////////////////////////////////
// Extracts number of CUDA cores from device property
{
    int cores = 0;
    int mp = devProp.multiProcessorCount;
    switch (devProp.major){
     case 2: // Fermi
      if (devProp.minor == 1) cores = mp * 48;
      else cores = mp * 32;
      break;
     case 3: // Kepler
      cores = mp * 192;
      break;
     case 5: // Maxwell
      cores = mp * 128;
      break;
     case 6: // Pascal
      if ((devProp.minor == 1) || (devProp.minor == 2)) cores = mp * 128;
      else if (devProp.minor == 0) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 7: // Volta and Turing
      if ((devProp.minor == 0) || (devProp.minor == 5)) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 8: // Ampere
      if (devProp.minor == 0) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     default:
      printf("Unknown device type\n");
      break;
      }
    return cores;
}

////////////////////////////////////////////////////////////////////////

__host__ void throwWARNING(const std::string messg)

////////////////////////////////////////////////////////////////////////
// Throws warning message 'messg'
{
   cout << endl << "WARNING: " << messg << endl << endl;
}

////////////////////////////////////////////////////////////////////////

__host__ void throwERROR(const std::string messg)

////////////////////////////////////////////////////////////////////////
// Throws error message 'messg'
{
   cout << endl << "ERROR: " << messg << endl << endl;
   cudaDeviceReset();
   exit(1);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
