//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// HopList operations
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains memory and I/O operations on arrays stored in the HopList structure
// PST, last update 06/19/2017

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "hoplist.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void HLfromComponents(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize, int *nbas, int *opid, int *mptr, prec_typ *coef, prec_typ *mat)

////////////////////////////////////////////////////////////////////////
// Packs a C++ HopList with data from a Fortran HopList

{
   H.nrk = nrk;
   H.ndof = ndof;
   H.nops = nops;
   H.msize = msize;
   H.nbas = nbas;
   H.opid = opid;
   H.mptr = mptr;
   H.coef = coef;
   H.mat = mat;
}

////////////////////////////////////////////////////////////////////////

__host__ void NewHopListonHost(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize)

////////////////////////////////////////////////////////////////////////
// Allocates HopList struct on host

{
   H.nrk = nrk;
   H.ndof = ndof;
   H.nops = nops;
   H.msize = msize;
   H.nbas = new int [ndof];
   H.opid = new int [nrk*ndof];
   H.mptr = new int [nops];
   H.coef = new prec_typ [nrk];
   H.mat = new prec_typ [msize];
}

////////////////////////////////////////////////////////////////////////

__host__ void FlushHopListonHost(HopList &H)

////////////////////////////////////////////////////////////////////////
// Deallocates HopList struct on host

{
   H.nrk = 0;
   H.ndof = 0;
   H.nops = 0;
   H.msize = 0;
   delete[] H.nbas;
   delete[] H.opid;
   delete[] H.mptr;
   delete[] H.coef;
   delete[] H.mat;
}


////////////////////////////////////////////////////////////////////////

__host__ void NewHopListonDevice(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize)

////////////////////////////////////////////////////////////////////////
// Allocates HopList struct on device

{
   H.nrk = nrk;
   H.ndof = ndof;
   H.nops = nops;
   H.msize = msize;
   cudaMalloc(&H.nbas,ndof*sizeof(int));
   cudaMalloc(&H.opid,nrk*ndof*sizeof(int));
   cudaMalloc(&H.mptr,nops*sizeof(int));
   cudaMalloc(&H.coef,nrk*sizeof(prec_typ));
   cudaMalloc(&H.mat,msize*sizeof(prec_typ));
}

////////////////////////////////////////////////////////////////////////

__host__ void FlushHopListonDevice(HopList &H)

////////////////////////////////////////////////////////////////////////
// Deallocates HopList struct on device

{
   H.nrk = 0;
   H.ndof = 0;
   H.nops = 0;
   H.msize = 0;
   cudaFree(H.nbas);
   cudaFree(H.opid);
   cudaFree(H.mptr);
   cudaFree(H.coef);
   cudaFree(H.mat);
}

////////////////////////////////////////////////////////////////////////

__host__ void SendHopList2GPU(HopList &v, HopList &w)

////////////////////////////////////////////////////////////////////////
// Copies HopList 'v' on host -> 'w' on device

{
   NewHopListonDevice(w,v.nrk,v.ndof,v.nops,v.msize);
   cudaMemcpy(w.nbas,v.nbas,v.ndof*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.opid,v.opid,v.nrk*v.ndof*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.mptr,v.mptr,v.nops*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.coef,v.coef,v.nrk*sizeof(prec_typ),cudaMemcpyHostToDevice);
   cudaMemcpy(w.mat,v.mat,v.msize*sizeof(prec_typ),cudaMemcpyHostToDevice);
}

////////////////////////////////////////////////////////////////////////

__host__ void GetHopListFromGPU(HopList &v, HopList &w)

////////////////////////////////////////////////////////////////////////
// Copies HopList 'v' on host <- 'w' on device
// 'v' should be allocated on host before calling

{
   cudaMemcpy(v.nbas,w.nbas,v.ndof*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.opid,w.opid,v.nrk*v.ndof*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.mptr,w.mptr,v.nops*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.coef,w.coef,v.nrk*sizeof(prec_typ),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.mat,w.mat,v.msize*sizeof(prec_typ),cudaMemcpyDeviceToHost);
}

////////////////////////////////////////////////////////////////////////

__host__ void PrintHopList(const HopList &H)

////////////////////////////////////////////////////////////////////////
// Prints out the CP-format Hamiltonian from the pointers and matrices
// stored in struct HopList (residing on host only)

{
   int i,j,k,l,m;

// Print the matrices for each term in H (lower triangle only)
   printf("\n");
   for (k=0;k<H.nrk;k++){  // Loop over terms in H
       printf(" C(%*d): %+E \n",5,k+1,H.coef[k]); // Print the coefficient

       for (j=0;j<H.ndof;j++){  // Loop over DOF
           i=0,l=1;
           while (i<H.nbas[j]){ // Prints the lower triangular matrix
               printf("%*d %*d",4,j+1,4,l);

               for (m=0;m<l;m++){
                   printf("  %+E ",H.mat[H.mptr[H.opid[k*H.ndof+j]-1]-1+i]);
                   i++;
                }
               l++;
               printf("\n");
           }
        printf("\n");
       }
   }
   printf("\n");
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Device-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__device__ prec_typ GetHopListelement(const HopList H, const int irk, const int idof, const int ielem)

////////////////////////////////////////////////////////////////////////
// Retrieves a HopList matrix element from device memory

{
   return H.mat[H.mptr[H.opid[irk*H.ndof+idof]-1]-1+ielem];
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
