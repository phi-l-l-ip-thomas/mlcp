//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CUDA ALS kernels
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains the CUDA kernels needed for Alternating Least Squares

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include "cusolverDn.h"
#include "global.h"
#include "utils.h"
#include "matrix.h"
#include "cpvector.h"
#include "modvecvec.h"
#include "linalg.h"
#include "als.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Constructors
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: init(const CPvec &Q_h, const int &rrk)

////////////////////////////////////////////////////////////////////////
// Constructor for mvpvars

{
   rk = rrk;
   nD = Q_h.nD;
   ndof = Q_h.ndof;
   nmax = *max_element(Q_h.nbas,Q_h.nbas+Q_h.ndof);
   nbas = new int [rrk];

// Inner product types
   mvp.initnew(nmax,0,rrk,nmax,ndof);
   psk.initnew(nmax,0,rrk,rrk,1);
   pst.initnew(nmax,0,rrk,rrk,1);
   bjk.initnew(rrk,1,rrk,nmax,1);

//   mvp.init(nmax,0,rrk,nmax,ndof);
//   psk.init(nmax,0,rrk,rrk,1);
//   psk.init(TILE_SIZE,0,rrk,rrk,1); // TILED
//   pst.init(nD,0,rrk,rrk,1);
//   pst.init(nmax,0,rrk,rrk,1); // PER MODE ("A")
//   bjk.init(rrk,1,rrk,nmax,1);

   for (int i=0;i<rrk;i++){
       nbas[i]=Q_h.nbas[i];
   }

// Kernel dimensions
   dB.x=dB.y=dB.z=1;
   dG.x=dG.y=dG.z=1;
   dB1.x=BLOCK_SIZE_1D;
   dB1.y=dB1.z=1;
   dG1.x=(rrk+BLOCK_SIZE_1D-1)/BLOCK_SIZE_1D;
   dG1.y=dG1.z=1;
   dB2.x=dB2.y=BLOCK_SIZE_2D;
   dB2.z=1;
   dG2.x=dG2.y=(rrk+BLOCK_SIZE_2D-1)/BLOCK_SIZE_2D;
   dG2.z=1;

// ALS arrays
   NewMatrixonDevice(bjk_d,rrk,nmax);
   NewMatrixonDevice(BB_d,rrk,rrk);

// CUSOLVER device arrays
   int LWORK;
   cusolverDnHandle_t als;
   cusolverDnCreate(&als);
#ifdef USE_DOUBLES
   cusolverDnDgetrf_bufferSize(als,rrk,rrk,&BB_d.mat[0],rrk,&LWORK);
#else
   cusolverDnSgetrf_bufferSize(als,rrk,rrk,&BB_d.mat[0],rrk,&LWORK);
#endif
   cusolverDnDestroy(als);
   cudaMalloc(&WORK,LWORK*sizeof(prec_typ));
   cudaMalloc(&INFO,sizeof(int));
   cudaMalloc(&IPV,rrk*sizeof(int));

// Normalization device arrays
   cudaMalloc(&tmp_d,rrk*sizeof(prec_typ));
   cudaMalloc(&sum_d,sizeof(prec_typ));
   cudaMalloc(&pen_d,sizeof(prec_typ));
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: newmvp(const CPvec &Q_h, const int &rrk)

////////////////////////////////////////////////////////////////////////
// Constructor for intwvars, matrix-vector products container

{
   which = 1;
   init(Q_h,rrk);
// CP-vectors on device
   NewCPveconDevice(F_d,Q_h.nbas,Q_h.ndof,rrk);
   NewCPveconDevice(G_d,Q_h.nbas,Q_h.ndof,rrk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: newgep(const CPvec &Q_h, const int &rrk, const int &nnbloc)

////////////////////////////////////////////////////////////////////////
// Constructor for intwvars, Gram-Schmidt container

{
   which = 2;
   init(Q_h,rrk);
   nbloc = nnbloc;
// CP-vectors on device
   NewCPveconDevice(F_d,Q_h.nbas,Q_h.ndof,rrk);
   NewCPveconDevice(G_d,Q_h.nbas,Q_h.ndof,rrk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: newgs(const CPvec &Q_h, const int &rrk, const int &nnbloc)

////////////////////////////////////////////////////////////////////////
// Constructor for intwvars, get Q^T H Q container

{
   which = 3;
   init(Q_h,rrk);
   nbloc = nnbloc;
   NewMatrixonDevice(BBk_d,rrk,rrk);
// Gram-Schmidt coefs
   cudaMalloc(&coef,nnbloc*sizeof(prec_typ));
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: newupdates(const CPvec &Q_h, const int &rrk, const int &nnbloc)

////////////////////////////////////////////////////////////////////////
// Constructor for intwvars, vector updates container

{
   which = 4;
   init(Q_h,rrk);
   nbloc = nnbloc;
// Vector to be optimized
   NewCPveconDevice(F_d,Q_h.nbas,Q_h.ndof,rrk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: newals(const CPvec &Q_h, const int &rrk, const int &nnbloc)

////////////////////////////////////////////////////////////////////////
// Constructor for intwvars, ordinary als container

{
   which = 5;
   init(Q_h,rrk);
   nbloc = nnbloc;
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: send(const CPvec &Q_h, const int &ivec)

////////////////////////////////////////////////////////////////////////
// Copies a CP-vector in the block from the host to the intwvars 
// variable "F_d" on the device

{
   SendCPvecRterms2GPU(Q_h,F_d,ivec*rk,rk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: point(CPvec &Q_d, const int &ivec)

////////////////////////////////////////////////////////////////////////
// Points the intwvars variable "F_d" on device to a vector in the block
// Q_d on the device

{
   PointCPvecRtermsonGPU(Q_d,Fo_d,ivec*rk,rk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: get(CPvec &Q_h, const int &ivec)

////////////////////////////////////////////////////////////////////////
// Copies a CP-vector in the block from the intwvars variable "F_d" on 
// the device to the block Q_h on the host

{
   GetCPvecRtermsfromGPU(Q_h,F_d,ivec*rk);
}

////////////////////////////////////////////////////////////////////////

__host__ void intwvars :: flush()

////////////////////////////////////////////////////////////////////////
// Disposes of intwvars struct

{
// Common to all types
   FlushMatrixonDevice(bjk_d);
   FlushMatrixonDevice(BB_d);
   cudaFree(WORK);
   cudaFree(INFO);
   cudaFree(IPV);
   cudaFree(tmp_d);
   cudaFree(sum_d);
   cudaFree(pen_d);
   delete[] nbas;

// Type specific
   if (which == 1){
      FlushCPveconDevice(F_d);
      FlushCPveconDevice(G_d);
   } else if (which == 2){
      FlushCPveconDevice(F_d);
      FlushCPveconDevice(G_d);
   } else if (which == 3){
      FlushMatrixonDevice(BBk_d);
      cudaFree(coef);
   } else if (which == 4){
      FlushCPveconDevice(F_d);
   } else if (which == 5){
      // Nothing to flush here
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


////////////////////////////////////////////////////////////////////////

__global__ void setpenalty_kernel(prec_typ *pen, const prec_typ *v, const int n)

////////////////////////////////////////////////////////////////////////
// Sets regularization penalty, device kernel

{

#ifdef USE_DOUBLES
   const prec_typ valpen = 1.0E-10;
#else
   const prec_typ valpen = 1.0E-5;
#endif

   int i;
   pen[0]=v[0];

   if (n > 1){
      for (i=1;i<n;i++){
          pen[0] = (v[i] > pen[0]) ? v[i] : pen[0];
      }
   }
   pen[0]*=valpen;
}

////////////////////////////////////////////////////////////////////////

__global__ void addpenalty_kernel(Matrix M, const prec_typ *E)

////////////////////////////////////////////////////////////////////////
// Adds ALS penalty to diagonal elements of matrix. No check is made to
// ensure that the matrix is square

{
// irk=index over diagonal elements of M
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

// Each thread adds 'E' to a diagonal element
   if (irk < M.nr && irk < M.nc){
      AddtoMatrixelement(M,irk,irk,E[0]);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void addpenalty_kernel(Matrix M, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Adds ALS penalty to diagonal elements of matrix. No check is made to
// ensure that the matrix is square

{
// irk=index over diagonal elements of M
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

// Each thread adds 'E' to a diagonal element
   if (irk < M.nr && irk < M.nc){
      AddtoMatrixelement(M,irk,irk,E);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void replacevwithsolns_kernel(CPvec v, const Matrix bjk, const int idof)

////////////////////////////////////////////////////////////////////////
// Replaces the base of a CP-vec along mode idof with solutions in bjk,
// device kernel

{
   int rk = v.nrk;
   int n = v.nbas[idof];

// irk=rank index, ibas=basis index of 'v', mode 'idof'
   int irk = blockIdx.x * blockDim.x + threadIdx.x;
   int ibas = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread multiplies a CP-element from v by 'E'
   if (irk < rk && ibas < n){
      SetCPelement(v,irk,idof,ibas,bjk.mat[irk+ibas*rk]);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k)

////////////////////////////////////////////////////////////////////////
// Constructs matrix bjk, the matrix of right-hand-sides of the linear
// system to be solved in ALS

{
   int j;
   int n = w.nbas[k];
   int rkv = bjk.nr;
   int rkw = w.nrk;
   prec_typ tmp;

// i1=row index, i2=col index, of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Add the next rkw terms to bjk
   if (i1 < rkv && i2 < n){
      tmp=GetMatrixelement(bjk,i1,i2); // Matrix element before adding the new terms
      for (j=0;j<rkw;j++){
          tmp+=w.coef[j]*w.base[j*w.nrdim[w.ndof]+w.nrdim[k]+i2]*PS.mat[j+i1*rkw];
      }
      SetMatrixelement(bjk,i1,i2,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const prec_typ *fac)

////////////////////////////////////////////////////////////////////////
// Constructs matrix bjk, the matrix of right-hand-sides of the linear
// system to be solved in ALS. This version multiplies the result by fac

{
   int j;
   int n = w.nbas[k];
   int rkv = bjk.nr;
   int rkw = w.nrk;
   prec_typ tmp;

// i1=row index, i2=col index, of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Add the next rkw terms to bjk
   if (i1 < rkv && i2 < n){
      for (j=0,tmp=0.0;j<rkw;j++){
          tmp+=w.coef[j]*w.base[j*w.nrdim[w.ndof]+w.nrdim[k]+i2]*PS.mat[j+i1*rkw];
      }
      tmp=fac[0]*tmp+GetMatrixelement(bjk,i1,i2); // Update matrix element
      SetMatrixelement(bjk,i1,i2,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_shared_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k)

////////////////////////////////////////////////////////////////////////
// Constructs matrix bjk, the matrix of right-hand-sides of the linear
// system to be solved in ALS

{
   int i3,j;
   int n = w.nbas[k];
   int rkv = bjk.nr;
   int rkw = w.nrk;
   prec_typ tmp;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[rkw];
   prec_typ* pm = (prec_typ*) &ws[rkw*blockDim.y];

// i1=row index, i2=col index, of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Load the coefficients of vector 'w' into shared memory
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;
   i3=cidx;
   while (i3 < rkw){
       wc[i3]=w.coef[i3];
       i3+=bloc;
   }
   __syncthreads();

// Load the base of the k-th dimension of vector 'w' into shared memory
   if (i2 < n){
      i3=threadIdx.x;
      while (i3 < rkw){
          tmp=GetCPelement(w,i3,k,i2);
          ws[i3+threadIdx.y*rkw]=tmp;
          i3+=blockDim.x;
      }
   }
   __syncthreads();

// Load columns of PS into shared memory
   if (i1 < rkv){
      i3=threadIdx.y;
      while (i3 < rkw){
          tmp=GetMatrixelement(PS,i3,i1);
          pm[i3+threadIdx.x*rkw]=tmp;
          i3+=blockDim.y;
      }
   }
   __syncthreads();

// Add the next rkw terms to bjk
   if (i1 < rkv && i2 < n){
      int s1 = threadIdx.x*rkw;
      int s2 = threadIdx.y*rkw;

      tmp=GetMatrixelement(bjk,i1,i2); // Matrix element before adding the new terms

      for (j=0;j<rkw;j++){
          tmp+=wc[j]*ws[j+s2]*pm[j+s1];
      }

      SetMatrixelement(bjk,i1,i2,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_shared_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const prec_typ *fac)

////////////////////////////////////////////////////////////////////////
// Constructs matrix bjk, the matrix of right-hand-sides of the linear
// system to be solved in ALS

{
   int i3,j;
   int n = w.nbas[k];
   int rkv = bjk.nr;
   int rkw = w.nrk;
   prec_typ tmp;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[rkw];
   prec_typ* pm = (prec_typ*) &ws[rkw*blockDim.y];

// i1=row index, i2=col index, of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Load the coefficients of vector 'w' into shared memory
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;
   i3=cidx;
   while (i3 < rkw){
       wc[i3]=w.coef[i3];
       i3+=bloc;
   }
   __syncthreads();

// Load the base of the k-th dimension of vector 'w' into shared memory
   if (i2 < n){
      i3=threadIdx.x;
      while (i3 < rkw){
          tmp=GetCPelement(w,i3,k,i2);
          ws[i3+threadIdx.y*rkw]=tmp;
          i3+=blockDim.x;
      }
   }
   __syncthreads();

// Load columns of PS into shared memory
   if (i1 < rkv){
      i3=threadIdx.y;
      while (i3 < rkw){
          tmp=GetMatrixelement(PS,i3,i1);
          pm[i3+threadIdx.x*rkw]=tmp;
          i3+=blockDim.y;
      }
   }
   __syncthreads();

// Add the next rkw terms to bjk
   if (i1 < rkv && i2 < n){
      int s1 = threadIdx.x*rkw;
      int s2 = threadIdx.y*rkw;

      for (j=0,tmp=0.0;j<rkw;j++){
          tmp+=wc[j]*ws[j+s2]*pm[j+s1];
      }
      tmp=fac[0]*tmp+GetMatrixelement(bjk,i1,i2); // Update matrix element
      SetMatrixelement(bjk,i1,i2,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_shared_kernelB(const CPvec w, Matrix PS, Matrix bjk, const int k, const int tilelen)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, mode k, shared device kernel with tiling.

{
   int tos;
   int j,i3,s1,s2;
   int ndof=w.ndof;
   int nrdim=w.nrdim[ndof];
   int vrk=bjk.nr;
   int wrk=w.nrk;
   int n=w.nbas[k];
   int ntile=(wrk+tilelen-1)/tilelen;
   int mos=w.nrdim[k]; // mode offset

   prec_typ acc;
   prec_typ tmp;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[tilelen];
   prec_typ* pm = (prec_typ*) &ws[tilelen*blockDim.y];

// i1=row index (v-rank), i2=col index (n), of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;
// All-threads index and count for loading coefs
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;

// Single-mode accumulator
   acc=GetMatrixelement(bjk,i1,i2);

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){

       tos=itile*tilelen; // tile offset

//     Load coefficients of vector 'w' into shared memory
       for (i3=cidx;i3<tilelen;i3+=bloc){
           tmp=0.0;
           j=tos+i3;
           if (j < wrk) tmp=w.coef[j];
	   wc[i3]=tmp;
       }
       __syncthreads();

//     Load base of w into shared memory
       for (i3=threadIdx.x;i3<tilelen;i3+=blockDim.x){
           tmp=0.0;
           j=tos+i3;
           if (i2 < n && j < wrk) tmp=w.base[j*nrdim+mos+i2];
           ws[i3+threadIdx.y*tilelen]=tmp;
       }
       __syncthreads();

//     Load tile of PS into shared memory
       for (i3=threadIdx.y;i3<tilelen;i3+=blockDim.y){
           tmp=0.0;
           j=tos+i3;
           if (i1 < vrk && j < wrk) tmp=PS.mat[j+i1*wrk];
           pm[i3+threadIdx.x*tilelen]=tmp;
       }
       __syncthreads();

//     Compute the inner product: loop over tile length
       for (i3=0,s1=threadIdx.x*tilelen,s2=threadIdx.y*tilelen;i3<tilelen;i3++,s1++,s2++){
           acc+=wc[i3]*pm[s1]*ws[s2];
       }
       __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   if (i1 < vrk && i2 < n) SetMatrixelement(bjk,i1,i2,acc);

}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_shared_kernelB(const CPvec w, Matrix PS, Matrix bjk, const int k, const int tilelen, const prec_typ *fac)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, mode k, shared device kernel with tiling.

{
   int tos;
   int j,i3,s1,s2;
   int ndof=w.ndof;
   int nrdim=w.nrdim[ndof];
   int vrk=bjk.nr;
   int wrk=w.nrk;
   int n=w.nbas[k];
   int ntile=(wrk+tilelen-1)/tilelen;
   int mos=w.nrdim[k]; // mode offset

   prec_typ acc;
   prec_typ tmp;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[tilelen];
   prec_typ* pm = (prec_typ*) &ws[tilelen*blockDim.y];

// i1=row index (v-rank), i2=col index (n), of bjk
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;
// All-threads index and count for loading coefs
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;

// Single-mode accumulator
   acc=0.0;

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){

       tos=itile*tilelen; // tile offset

//     Load coefficients of vector 'w' into shared memory
       for (i3=cidx;i3<tilelen;i3+=bloc){
           tmp=0.0;
           j=tos+i3;
           if (j < wrk) tmp=w.coef[j];
           wc[i3]=tmp;
       }
       __syncthreads();

//     Load base of w into shared memory
       for (i3=threadIdx.x;i3<tilelen;i3+=blockDim.x){
           tmp=0.0;
           j=tos+i3;
           if (i2 < n && j < wrk) tmp=w.base[j*nrdim+mos+i2];
           ws[i3+threadIdx.y*tilelen]=tmp;
       }
       __syncthreads();

//     Load tile of PS into shared memory
       for (i3=threadIdx.y;i3<tilelen;i3+=blockDim.y){
           tmp=0.0;
           j=tos+i3;
           if (i1 < vrk && j < wrk) tmp=PS.mat[j+i1*wrk];
           pm[i3+threadIdx.x*tilelen]=tmp;
       }
       __syncthreads();

//     Compute the inner product: loop over tile length
       for (i3=0,s1=threadIdx.x*tilelen,s2=threadIdx.y*tilelen;i3<tilelen;i3++,s1++,s2++){
           acc+=wc[i3]*pm[s1]*ws[s2];
       }
       __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   if (i1 < vrk && i2 < n) {
      tmp=fac[0]*acc+GetMatrixelement(bjk,i1,i2); // Update matrix element
      SetMatrixelement(bjk,i1,i2,tmp);
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_registerized_kernel(const CPvec w, Matrix PS, Matrix bjk, const int idof, const int tilelen)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different terms in
// the CP expansion, mode idof, shared device kernel with registerized loads.

{

   volatile int ilg,irg,icg,irow,icol,irk,iin;
   volatile int tos;
   int i,j,k;
   int ndof=w.ndof;
   int nrdim=w.nrdim[ndof];
   int vrk=bjk.nr;
   int wrk=w.nrk;
   int n=w.nbas[idof];
   int ntile=(wrk+tilelen-1)/tilelen;
   int mos=w.nrdim[idof]; // mode offset
   prec_typ acc[WPT_X*WPT_Y];
   prec_typ wreg[WPT_Y];
   prec_typ creg;
   prec_typ preg;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[tilelen];
   prec_typ* pm = (prec_typ*) &ws[tilelen*blockDim.y*WPT_Y];

// All-threads index and count for loading coefs
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;
// Block offsets ({x,y} for {rank-of-v,n}, respectively)
   int bosx = blockIdx.x * blockDim.x * WPT_X;
   int bosy = blockIdx.y * blockDim.y * WPT_Y;
// column-load-groups
   int clgx = (tilelen + blockDim.y - 1) / blockDim.y;
   int clgy = (tilelen + blockDim.x - 1) / blockDim.x;
// row-sets
   int rsx = (tilelen + blockDim.y - 1) / tilelen;
   int rsy = (tilelen + blockDim.x - 1) / tilelen;
// row-load-groups
   int rlgx = (WPT_X + rsx - 1) / rsx;
   int rlgy = (WPT_Y + rsy - 1) / rsy;
// total number of load-groups
   int nlgx = clgx*rlgx;
   int nlgy = clgy*rlgy;
// index for row-set
   int irsx = threadIdx.y / tilelen;
   int irsy = threadIdx.x / tilelen;

// Single-mode product registers
   for (i=0;i<WPT_X*WPT_Y;i++) acc[i]=0.0;

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){

       tos=itile*tilelen; // tile offset

//     Load coefficients of vector 'w' into shared memory
       for (i=cidx;i<tilelen;i+=bloc){
           j=tos+i;
           creg=0.0;
           if (j < wrk) creg=w.coef[j];
           wc[i]=creg;
       }
       __syncthreads();

//    Load tile of w into shared memory
      for (ilg=0;ilg<nlgy;ilg++){
          irg=ilg/clgy; // row-group index (basis n)
          icg=ilg%clgy; // col-group index (rank-of-w)
          irow=threadIdx.y+(irg+irsy*rlgy)*blockDim.y;
          icol=(threadIdx.x%tilelen)+icg*blockDim.x;
          iin=irow+bosy; // w basis
	  irk=tos+icol;  // w rank
          creg=0.0;
          if (irk < wrk && iin < n) creg=w.base[irk*nrdim+mos+iin];
          ws[icol+irow*tilelen]=creg;
      }
      __syncthreads();

//     Load tile of PS into shared memory
       for (ilg=0;ilg<nlgx;ilg++){
           irg=ilg/clgx; // row-group index (rank-of-v)
           icg=ilg%clgx; // col-group index (rank-of-w)
           irow=threadIdx.x+(irg+irsx*rlgx)*blockDim.x;
           icol=(threadIdx.y%tilelen)+icg*blockDim.y;
           irk=irow+bosx; // v-rank in PS
	   iin=tos+icol;  // w-rank in PS
           creg=0.0;
	   if (irk < vrk && iin < wrk) creg=PS.mat[iin+irk*wrk];
	   pm[icol+irow*tilelen]=creg;
       }
       __syncthreads();

//    Compute the inner product: loop over tile length
      for (i=0;i<tilelen;i++){

          creg=wc[i]; // Coef of w

//        Cache entries of w into registers and premultiply coef
          for (k=0;k<WPT_Y;k++){
              icol=threadIdx.y+k*blockDim.y;
              wreg[k]=creg*ws[i+icol*tilelen];
          }

//        Compute the inner product
          for (j=0;j<WPT_X;j++){
              irow=threadIdx.x+j*blockDim.x;
              preg=pm[i+irow*tilelen]; // Cached PS entry
              for (k=0;k<WPT_Y;k++){
                  acc[j+WPT_X*k]+=preg*wreg[k];
              }
          }

      }
      __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   for (j=0;j<WPT_X;j++){
       irow=bosx+threadIdx.x+j*blockDim.x;
       for (k=0;k<WPT_Y;k++){
           icol=bosy+threadIdx.y+k*blockDim.y;
	   i=j+WPT_X*k;
           if (irow < vrk && icol < n) {
              creg=acc[i]+GetMatrixelement(bjk,irow,icol);
              SetMatrixelement(bjk,irow,icol,creg);
           }

       }
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constbjk_registerized_kernel(const CPvec w, Matrix PS, Matrix bjk, const int idof, const int tilelen, const prec_typ *fac)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different terms in
// the CP expansion, mode idof, shared device kernel with registerized loads.

{

   volatile int ilg,irg,icg,irow,icol,irk,iin;
   volatile int tos;
   int i,j,k;
   int ndof=w.ndof;
   int nrdim=w.nrdim[ndof];
   int vrk=bjk.nr;
   int wrk=w.nrk;
   int n=w.nbas[idof];
   int ntile=(wrk+tilelen-1)/tilelen;
   int mos=w.nrdim[idof]; // mode offset
   prec_typ acc[WPT_X*WPT_Y];
   prec_typ wreg[WPT_Y];
   prec_typ creg;
   prec_typ preg;

   extern __shared__ prec_typ wscp[];
   prec_typ* wc = (prec_typ*) wscp;
   prec_typ* ws = (prec_typ*) &wc[tilelen];
   prec_typ* pm = (prec_typ*) &ws[tilelen*blockDim.y*WPT_Y];

// All-threads index and count for loading coefs
   int cidx = blockDim.x * threadIdx.y + threadIdx.x;
   int bloc = blockDim.x * blockDim.y;
// Block offsets ({x,y} for {rank-of-v,n}, respectively)
   int bosx = blockIdx.x * blockDim.x * WPT_X;
   int bosy = blockIdx.y * blockDim.y * WPT_Y;
// column-load-groups
   int clgx = (tilelen + blockDim.y - 1) / blockDim.y;
   int clgy = (tilelen + blockDim.x - 1) / blockDim.x;
// row-sets
   int rsx = (tilelen + blockDim.y - 1) / tilelen;
   int rsy = (tilelen + blockDim.x - 1) / tilelen;
// row-load-groups
   int rlgx = (WPT_X + rsx - 1) / rsx;
   int rlgy = (WPT_Y + rsy - 1) / rsy;
// total number of load-groups
   int nlgx = clgx*rlgx;
   int nlgy = clgy*rlgy;
// index for row-set
   int irsx = threadIdx.y / tilelen;
   int irsy = threadIdx.x / tilelen;

// Single-mode product registers
   for (i=0;i<WPT_X*WPT_Y;i++) acc[i]=0.0;

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){

       tos=itile*tilelen; // tile offset

//     Load coefficients of vector 'w' into shared memory
       for (i=cidx;i<tilelen;i+=bloc){
           j=tos+i;
           creg=0.0;
           if (j < wrk) creg=w.coef[j];
           wc[i]=creg;
       }
       __syncthreads();

//    Load tile of w into shared memory
      for (ilg=0;ilg<nlgy;ilg++){
          irg=ilg/clgy; // row-group index (basis n)
          icg=ilg%clgy; // col-group index (rank-of-w)
          irow=threadIdx.y+(irg+irsy*rlgy)*blockDim.y;
          icol=(threadIdx.x%tilelen)+icg*blockDim.x;
          iin=irow+bosy; // w basis
	  irk=tos+icol;  // w rank
          creg=0.0;
          if (irk < wrk && iin < n) creg=w.base[irk*nrdim+mos+iin];
          ws[icol+irow*tilelen]=creg;
      }
      __syncthreads();

//     Load tile of PS into shared memory
       for (ilg=0;ilg<nlgx;ilg++){
           irg=ilg/clgx; // row-group index (rank-of-v)
           icg=ilg%clgx; // col-group index (rank-of-w)
           irow=threadIdx.x+(irg+irsx*rlgx)*blockDim.x;
           icol=(threadIdx.y%tilelen)+icg*blockDim.y;
           irk=irow+bosx; // v-rank in PS
	   iin=tos+icol;  // w-rank in PS
           creg=0.0;
	   if (irk < vrk && iin < wrk) creg=PS.mat[iin+irk*wrk];
	   pm[icol+irow*tilelen]=creg;
       }
       __syncthreads();

//    Compute the inner product: loop over tile length
      for (i=0;i<tilelen;i++){

          creg=wc[i]; // Coef of w

//        Cache entries of w into registers and premultiply coef
          for (k=0;k<WPT_Y;k++){
              icol=threadIdx.y+k*blockDim.y;
              wreg[k]=creg*ws[i+icol*tilelen];
          }

//        Compute the inner product
          for (j=0;j<WPT_X;j++){
              irow=threadIdx.x+j*blockDim.x;
              preg=pm[i+irow*tilelen]; // Cached PS entry
              for (k=0;k<WPT_Y;k++){
                  acc[j+WPT_X*k]+=preg*wreg[k];
              }
          }

      }
      __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   for (j=0;j<WPT_X;j++){
       irow=bosx+threadIdx.x+j*blockDim.x;
       for (k=0;k<WPT_Y;k++){
           icol=bosy+threadIdx.y+k*blockDim.y;
	   i=j+WPT_X*k;
	   if (irow < vrk && icol < n) {
              creg=fac[0]*acc[i]+GetMatrixelement(bjk,irow,icol);
              SetMatrixelement(bjk,irow,icol,creg);
           }
       }
   }

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void als_build_linear_system(const PVVkernel &pst, const PVVkernel &bjk, CPvec &F_d, CPvec &Q_d, Matrix &BB_d, Matrix &bjk_d, const int &idof, const int &ni, const int &iv)

////////////////////////////////////////////////////////////////////////
// Builds the BB, bjk matrices required for solving the linear systems
// for the intertwined power method.

{
   int k;
   CPvec G_d;

// Initialize values
   int rk = F_d.nrk;
   prec_typ beta = 0.0;

// Initialize CUDA types
   dim3 dB1(BLOCK_SIZE_1D);
   dim3 dB2(BLOCK_SIZE_2D,BLOCK_SIZE_2D);
   dim3 dG1((rk+dB1.x-1)/dB1.x);                           // rk x 1
   dim3 dimGridbjk((rk+dB2.x-1)/dB2.x,(ni+dB2.y-1)/dB2.y); // rk x n

// Initialize bjk on device
   fillmatrix_kernel<<<dimGridbjk,dB2>>>(bjk_d,beta);

// Build the r.h.s. of the linear system. Here we loop over rk-sized 
// blocks of terms in G. The k index begins at 1 (not 0) since the first
// rk terms belong to F
   for (k=1;k<iv;k++){

//     Designate the k-th block of terms in Q as "G_d"
       PointCPvecRtermsonGPU(Q_d,G_d,k*rk,rk);

//     Construct PS(!=idof) matrix
       constpstot_wrapper(pst,G_d,F_d,BB_d,idof,false,false);

//     Add terms in PS(!=idof) to bjk matrix
       constbjk_wrapper(bjk,G_d,BB_d,bjk_d,idof);
   }

// Build the l.h.s. of the linear system, BB(!=idof)
   constpstot_wrapper(pst,F_d,F_d,BB_d,idof,false,false);

}

////////////////////////////////////////////////////////////////////////

__host__ void ITERATE_ALS(CPvec &Q_h, CPvec &Q_d, const int &rk, const int &nbloc, const int &nals)

////////////////////////////////////////////////////////////////////////
// Vanilla ALS driver. Input Q_d contains F as first rk terms and G as
// the remaining rk*(nbloc-1) terms

{
   if (nals == 0) return;

// Initialize ALS container
   intwvars IV;
   IV.newals(Q_h,rk,nbloc);
   CPvec F_d;
   int ndof = Q_d.ndof;
   int nD = IV.nD;

// Initialize CUDA types
   dim3 dimGridb((nbloc+IV.dB1.x-1)/IV.dB1.x); // nbloc x 1
   dim3 dimGridCopy((rk+IV.dB2.x-1)/IV.dB2.x,(nD+IV.dB2.y-1)/IV.dB2.y); // rk x nD

// Designate the first rk terms in Q as "F_d"
   PointCPvecRtermsonGPU(Q_d,F_d,0,rk);

// Compute the regularization penalty from the coefs of F
   setpenalty_kernel<<<IV.dG,IV.dB>>>(IV.pen_d,F_d.coef,rk);

// Loop over intertwining interations
   for (int i=0;i<nals;i++){

//     Loop over DOF
       for (int j=0;j<ndof;j++){
           int n = IV.nbas[j];

//         Calculate BB and b_j_k
	   als_build_linear_system(IV.pst,IV.bjk,F_d,Q_d,IV.BB_d,IV.bjk_d,j,n,nbloc);

//         Add penalty to BB to avoid ill-conditioning
           addpenalty_kernel<<<IV.dG1,IV.dB1>>>(IV.BB_d,IV.pen_d);

//         Solve linear system B*c_j_k = b_j_k
           solve_linsys(IV.BB_d,IV.bjk_d,IV.WORK,IV.IPV,IV.INFO,rk,n);

//         Construct improved F and normalize base
           dim3 dimGridbjk((rk+IV.dB2.x-1)/IV.dB2.x,(n+IV.dB2.y-1)/IV.dB2.y); // rF x n
           replacevwithsolns_kernel<<<dimGridbjk,IV.dB2>>>(F_d,IV.bjk_d,j);
           normonebase_kernel<<<IV.dG1,IV.dB1>>>(F_d,j);

       } // Loop over DOF j

   } // Loop over iterations i

// Copy the first rk terms of Q_d (the ones corresponding to F) to host
   GetCPvecRtermsfromGPU(Q_h,F_d,0);

// Dispose ALS container
   IV.flush();

}

////////////////////////////////////////////////////////////////////////

__host__ void constbjk_wrapper(const PVVkernel &PV, const CPvec &w, const Matrix &PS, Matrix &bjk, const int &k)

////////////////////////////////////////////////////////////////////////
// Wrapper for shared/unshared calls to constbjk kernel

{
   if (PV.shared == 2) {
      constbjk_registerized_kernel<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k,PV.tilelen);
   } else if (PV.shared == 1) {
//      constbjk_shared_kernel<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k);
      constbjk_shared_kernelB<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k,PV.tilelen);
   } else {
      constbjk_kernel<<<PV.dG,PV.dB>>>(w,PS,bjk,k);
   }
}

////////////////////////////////////////////////////////////////////////

__host__ void constbjk_wrapper(const PVVkernel &PV, const CPvec &w, const Matrix &PS, Matrix &bjk, const int &k, const prec_typ *fac)

////////////////////////////////////////////////////////////////////////
// Wrapper for shared/unshared calls to constbjk kernel
// This version multiplies the result by 'fac'

{
   if (PV.shared == 2) {
      constbjk_registerized_kernel<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k,PV.tilelen,fac);
   } else if (PV.shared == 1) {
      constbjk_shared_kernel<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k,fac);
//      constbjk_shared_kernelB<<<PV.dG,PV.dB,PV.mem>>>(w,PS,bjk,k,PV.tilelen,fac);
   } else {
      constbjk_kernel<<<PV.dG,PV.dB>>>(w,PS,bjk,k,fac);
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
