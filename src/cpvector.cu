//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CP-vector operations
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains memory and I/O operations on arrays stored in CP data-compressed format
// and belonging to class CPvec
// PST, last update 05/30/2017

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "cpvector.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__global__ void copycpcoef_kernel(const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Copies a CP-vec coef from one CP-vec to another, device kernel
// v and w must be the same size, as no check is performed

{
// irk=rank index of v
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk) w.coef[irk] = v.coef[irk];
}

////////////////////////////////////////////////////////////////////////

__global__ void copycpbase_kernel(const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Copies a CP-vec base from one CP-vec to another, device kernel
// v and w must be the same size, as no check is performed

{
   int rk = v.nrk;
   int nD = v.nD;

// irk=rank index, ibD=base index (covers all D DOFs)
   int irk = blockIdx.x * blockDim.x + threadIdx.x;
   int ibD = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread extracts a CP-element from v and places in w
   if (irk < rk && ibD < nD){
      int j = irk*nD+ibD;
      w.base[j] = v.base[j];
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void setcpcoef_kernel(CPvec v, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Assigns a CP-vec coef to value E, device kernel

{
// irk=rank index of v
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk) v.coef[irk] = E;
}

////////////////////////////////////////////////////////////////////////

__global__ void scalecpcoef_kernel(CPvec v, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Multiplies a CP-vec coef in device memory by a POSITIVE value. Since
// coefs are always positive, one must also call changecpsign_kernel if
// multiplication by a negative value is desired

{
// irk=rank index of v
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk) v.coef[irk] *= abs(E);
}

////////////////////////////////////////////////////////////////////////

__global__ void scalecpbase_kernel(CPvec v, const int idof, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Multiplies the base of a CP-vec by 'E' along mode idof, device kernel

{
   int rk = v.nrk;
   int n = v.nbas[idof];

// irk=rank index, ibas=basis index of 'v', mode 'idof'
   int irk = blockIdx.x * blockDim.x + threadIdx.x;
   int ibas = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread multiplies a CP-element from v by 'E'
   if (irk < rk && ibas < n){
      MultiplyCPelement(v,irk,idof,ibas,E);
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void CPfromComponents(CPvec &v, prec_typ *base, prec_typ *coef, int *nbas, const int &ndof, const int &nrk)

////////////////////////////////////////////////////////////////////////
// Wrapper for initializing a CP-vector using the C++ class using data
// from a Fortran CP-vector
{
   int i;
   int nrdim;

   v.nrk=nrk;
   v.ndof=ndof;
   v.base=base;
   v.coef=coef;
   v.nbas=nbas;
   v.nmax=0;

   v.nrdim = new int [ndof+1];
   for (i=0,nrdim=0;i<ndof;i++){
       if (nbas[i]<1) throw("nbas must be positive");
       if (nbas[i]>v.nmax) v.nmax=nbas[i];
       v.nrdim[i]=nrdim;
       nrdim=nrdim+nbas[i];
   }
   v.nrdim[ndof]=v.nD=nrdim;

}

////////////////////////////////////////////////////////////////////////

__host__ void NewCPveconHost(CPvec &v, const int *nbas, const int &ndof, const int &nrk)

////////////////////////////////////////////////////////////////////////
// Initializes CP-vector on host

{
   int i;
   int nrdim;

   v.nrk = nrk;
   v.ndof = ndof;
   v.nbas = new int [ndof];
   v.nrdim = new int [ndof+1];
   v.coef = new prec_typ [nrk];
   v.nmax = 0;

   for (i=0,nrdim=0;i<ndof;i++){
       if (nbas[i]<1) throw("nbas must be positive");
       if (nbas[i]>v.nmax) v.nmax=nbas[i];
       v.nrdim[i]=nrdim;
       nrdim=nrdim+nbas[i];
       v.nbas[i]=nbas[i];
   }

// Define v.nrdim[ndof] AND v.nD since it is useful to have this value
// as an array element in device memory and as a parameter on the host
   v.nrdim[ndof]=v.nD=nrdim;
   v.base = new prec_typ [nrdim*nrk];

}

////////////////////////////////////////////////////////////////////////

__host__ void FlushCPveconHost(CPvec &v)

////////////////////////////////////////////////////////////////////////
// Deallocates CP-vector on host

{
   v.nrk = 0;
   v.ndof = 0;
   v.nD = 0;
   delete[] v.nbas;
   delete[] v.nrdim;
   delete[] v.coef;
   delete[] v.base;

}

////////////////////////////////////////////////////////////////////////

__host__ void NewCPveconDevice(CPvec &v, const int *nbas, const int &ndof, const int &nrk)

////////////////////////////////////////////////////////////////////////
// Initializes/allocates CP-vector on device

{
   int i;
   int nrdm;
   int *nrdim;

   v.nrk = nrk;
   v.ndof = ndof;
   v.nmax = 0;
   nrdim = new int [ndof+1];
   for (i=0,nrdm=0;i<ndof;i++){
       if (nbas[i]<1) throw("nbas must be positive");
       if (nbas[i]>v.nmax) v.nmax=nbas[i];
       nrdim[i]=nrdm;
       nrdm+=nbas[i];
   }
   nrdim[ndof] = nrdm;
   v.nD=nrdm;

   cudaMalloc(&v.nbas,ndof*sizeof(int));
   cudaMemcpy(v.nbas,nbas,ndof*sizeof(int),cudaMemcpyHostToDevice);
   cudaMalloc(&v.nrdim,(ndof+1)*sizeof(int));
   cudaMemcpy(v.nrdim,nrdim,(ndof+1)*sizeof(int),cudaMemcpyHostToDevice);
   cudaMalloc(&v.coef,nrk*sizeof(prec_typ));
   cudaMalloc(&v.base,nrdm*nrk*sizeof(prec_typ));

   delete[] nrdim;
}

////////////////////////////////////////////////////////////////////////

__host__ void FlushCPveconDevice(CPvec &v)

////////////////////////////////////////////////////////////////////////
// Deallocates CP-vector on device

{
   v.nrk = 0;
   v.ndof = 0;
   v.nD = 0;
   v.nmax = 0;
   cudaFree(v.nbas);
   cudaFree(v.nrdim);
   cudaFree(v.coef);
   cudaFree(v.base);
}

////////////////////////////////////////////////////////////////////////

__host__ void SendCPvec2GPU(const CPvec &v, CPvec &w)

////////////////////////////////////////////////////////////////////////
// Copies CP-vector 'v' on host -> 'w' on device

{
   NewCPveconDevice(w,v.nbas,v.ndof,v.nrk);
   cudaMemcpy(w.nbas,v.nbas,v.ndof*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.nrdim,v.nrdim,(v.ndof+1)*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.coef,v.coef,v.nrk*sizeof(prec_typ),cudaMemcpyHostToDevice);
   cudaMemcpy(w.base,v.base,v.nrdim[v.ndof]*v.nrk*sizeof(prec_typ),cudaMemcpyHostToDevice);
}

////////////////////////////////////////////////////////////////////////

__host__ void SendCPvecRterms2GPU(const CPvec &v, CPvec &w, const int &iR, const int &R)

////////////////////////////////////////////////////////////////////////
// Copies R terms of CP-vector 'v' on host -> 'w' on device, where 'R'
// is the number of terms in 'w', and 'iR' is the first term to copy

{
   if (iR < 0 || iR > v.nrk){
       cout << "iR = " << iR << endl;
       throw("SendCPvecRterms2GPU(): invalid initial rank index iR");
   }
   if (iR+R > v.nrk){
       cout << "iR+R = " << iR+R << "; v.nrk = " << v.nrk << endl;
       throw("SendCPvecRterms2GPU(): max rank to copy exceeds bounds");
   }

   int ibas=iR*v.nrdim[v.ndof];

   cudaMemcpy(w.nbas,v.nbas,v.ndof*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.nrdim,v.nrdim,(v.ndof+1)*sizeof(int),cudaMemcpyHostToDevice);
   cudaMemcpy(w.coef,&v.coef[iR],R*sizeof(prec_typ),cudaMemcpyHostToDevice);
   cudaMemcpy(w.base,&v.base[ibas],R*v.nrdim[v.ndof]*sizeof(prec_typ),cudaMemcpyHostToDevice);
}

////////////////////////////////////////////////////////////////////////

__host__ void PointCPvecRtermsonGPU(CPvec &v, CPvec &w, const int &iR, const int &R)

////////////////////////////////////////////////////////////////////////
// Sets 'w' on device to point to R terms of CP-vector 'v' on device,
// where 'R' is the number of terms in 'w', and 'iR' is the offset

{
   if (iR < 0 || iR > v.nrk){
       cout << "iR = " << iR << endl;
       throw("PointCPvecRtermsonGPU(): invalid initial rank index iR");
   }
   if (iR+R > v.nrk){
       cout << "iR+R = " << iR+R << "; v.nrk = " << v.nrk << endl;
       throw("PointCPvecRtermsonGPU(): max rank to copy exceeds bounds");
   }

   w.nrk=R;
   w.ndof=v.ndof;
   w.nbas=v.nbas;
   w.nrdim=v.nrdim;
   w.nD=v.nD;
   w.nmax=v.nmax;
   w.coef=&v.coef[iR];
   w.base=&v.base[iR*v.nD];
}

////////////////////////////////////////////////////////////////////////

__host__ void GetCPvecFromGPU(CPvec &v, const CPvec &w)

////////////////////////////////////////////////////////////////////////
// Copies CP-vector 'v' on host <- 'w' on device
// 'v' should be allocated on host before calling

{
   cudaMemcpy(v.nbas,w.nbas,v.ndof*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.nrdim,w.nrdim,(v.ndof+1)*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.coef,w.coef,v.nrk*sizeof(prec_typ),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.base,w.base,v.nrdim[v.ndof]*v.nrk*sizeof(prec_typ),cudaMemcpyDeviceToHost);
}

////////////////////////////////////////////////////////////////////////

__host__ void GetCPvecRtermsfromGPU(CPvec &v, const CPvec &w, const int &iR)

////////////////////////////////////////////////////////////////////////
// Overwrites R terms of CP-vector 'v' on host <- 'w' on device, where 
// 'iR' is the term in 'v' where overwriting begins

{
   if (iR < 0 || iR > v.nrk){
       cout << "iR = " << iR << endl;
       throw("GetCPvecRtermsfromGPU(): invalid initial rank index iR");
   }
   if (iR+w.nrk > v.nrk){
       cout << "iR+R = " << iR+w.nrk << "; v.nrk = " << v.nrk << endl;
       throw("GetCPvecRtermsfromGPU(): max rank to copy exceeds bounds");
   }

   int ibas=iR*v.nrdim[v.ndof];

   cudaMemcpy(v.nbas,w.nbas,v.ndof*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(v.nrdim,w.nrdim,(v.ndof+1)*sizeof(int),cudaMemcpyDeviceToHost);
   cudaMemcpy(&v.coef[iR],w.coef,w.nrk*sizeof(prec_typ),cudaMemcpyDeviceToHost);
   cudaMemcpy(&v.base[ibas],w.base,v.nrdim[v.ndof]*w.nrk*sizeof(prec_typ),cudaMemcpyDeviceToHost);

}

////////////////////////////////////////////////////////////////////////

__host__ void PrintCPvec(const CPvec &v)

////////////////////////////////////////////////////////////////////////
// Prints out CP-vector residing on host

{
   int i,j,k;

// Print the coefficient
   cout << endl << " Vcoef =  ";
   for (j=0;j<v.nrk;j++){
       printf("  %E ",v.coef[j]);
   }
   cout << endl;

// Print the base
   for (i=0;i<v.ndof;i++){
       for (k=0;k<v.nbas[i];k++){
           printf("%*d %*d ",4,i+1,4,k+1);
           for (j=0;j<v.nrk;j++){
               printf(" %13.10f ",v.base[j*v.nrdim[v.ndof]+v.nrdim[i]+k]);
           }
           cout << endl;
       }
   }
   cout << endl;
}

////////////////////////////////////////////////////////////////////////

__host__ int getsharedmemCP(const dim3 &dB, const int &n, const int &m)

////////////////////////////////////////////////////////////////////////
// Computes the amount of shared memory required by the shared version
// of kernels using CP-vecs

{
   return sizeof(prec_typ)*n*(dB.x+dB.y+m);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Device-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__device__ prec_typ GetCPelement(const CPvec v, const int irk, const int idof, const int ibas)

////////////////////////////////////////////////////////////////////////
// Retrieves a CP-vec element from device memory

{
   return v.base[irk*v.nrdim[v.ndof] + v.nrdim[idof] + ibas];
}

////////////////////////////////////////////////////////////////////////

__device__ void SetCPelement(CPvec v, const int irk, const int idof, const int ibas, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Assigns a CP-vec element in device memory

{
   v.base[irk*v.nrdim[v.ndof] + v.nrdim[idof] + ibas]=E;
}

////////////////////////////////////////////////////////////////////////

__device__ void MultiplyCPelement(CPvec v, const int irk, const int idof, const int ibas, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Multiplies a CP-vec element in device memory

{
   v.base[irk*v.nrdim[v.ndof] + v.nrdim[idof] + ibas]*=E;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
