//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CP-vector sums and products
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains subroutines for computing inner products of arrays stored in CP tensor
// data compressed format. Normalization routines are also found here.
// PST, last update 05/30/2017

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "matrix.h"
#include "cpvector.h"
#include "utils.h"
#include "modvecvec.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Constructors
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void PVVkernel :: initnew(const int &nn, const int &mm, const int &rrkx, const int &rrky, const int &rrkz)

////////////////////////////////////////////////////////////////////////
// Initializes PVVkernel

{
// Initialize values
   n = nn;
   m = mm;
   rkx = rrkx;
   rky = rrky;
   rkz = rrkz;
   shared = 1;

// Calculate the 'base size' of the block: contract block sizes by half until they are the
// smallest power of two that is greater than or equal to the ranks. This avoids scheduling 
// an excessive number of threads if the total block size is small
   int basex = BLOCK_SIZE_2D;
   while (basex >= 2*rkx){
       basex /= 2;
   }
   int basey = BLOCK_SIZE_2D;
   while (basey >= 2*rky){
       basey /= 2;
   }

// Make sure base sizes give a total block size of at least WARP_SIZE
   while (basex*basey < WARP_SIZE){
       basey *= 2;
       if (basex*basey >= WARP_SIZE) break;
       basex *= 2;
   }

// Set the block dimensions
   dB.x = basex;
   dB.y = basey;
   dB.z = 1;

// If the grid dimensions exceed the work-per-thread in BOTH X and Y dimensions,
// select the registerized code
   if ((rkx+dB.x-1)/dB.x >= WPT_X && (rky+dB.y-1)/dB.y >= WPT_Y){
//    Temporarily scale block dimensions by work-per-thread to compute memory
      dB.x = basex*WPT_X;
      dB.y = basey*WPT_Y;
      shared = 2;
   }

// Adjust the tile size so block sizes do not exceed shared memory
   tilelen = 1;
   while (getsharedmemCP(dB,2*tilelen,mm) <= SHARED_MEMORY){
       if (tilelen >= nn) break;
       tilelen *= 2;
   }

// Set the memory and grid dimensions
   mem = 0;
   if (shared > 0) mem = getsharedmemCP(dB,tilelen,mm);
   dG.x = (rkx+dB.x-1)/dB.x;  // x-dim over rank of 1st vector
   dG.y = (rky+dB.y-1)/dB.y;  // y-dim over rank of 2nd vector
   dG.z = (rkz+dB.z-1)/dB.z;  // z-dim over degrees-of-freedom

// Unscale block dimensions for correct launch parameters
   if (shared == 2){
      dB.x /= WPT_X;
      dB.y /= WPT_Y;
   }

}

////////////////////////////////////////////////////////////////////////

__host__ void PVVkernel :: init(const int &nn, const int &mm, const int &rrkx, const int &rrky, const int &rrkz)

////////////////////////////////////////////////////////////////////////
// Initializes PVVkernel

{
// Initialize values
   n = nn;
   m = mm;
   rkx = rrkx;
   rky = rrky;
   rkz = rrkz;
   shared = 1;

// Calculate the 'base size' of the block: contract block sizes by half until they are the
// smallest power of two that is greater than or equal to the ranks. This avoids scheduling 
// an excessive number of threads if the total block size is small
   int basex = BLOCK_SIZE_2D;
   while (basex >= 2*rkx){
       basex /= 2;
   }
   int basey = BLOCK_SIZE_2D;
   while (basey >= 2*rky){
       basey /= 2;
   }

// Make sure base sizes give a total block size of at least 32
   while (basex*basey < 32){
       basey *= 2;
       if (basex*basey >= 32) break;
       basex *= 2;
   }

// Set the block dimensions
   dB.x = basex;
   dB.y = basey;
   dB.z = 1;

// Make sure the block dimensions are small enough for the device shared memory
   while (getsharedmemCP(dB,nn,mm) > SHARED_MEMORY && dB.x*dB.y >= 32) {
      dB.y /= 2;
      if (getsharedmemCP(dB,nn,mm) <= SHARED_MEMORY || dB.x*dB.y < 32) break;
      dB.x /= 2;
   }
// If the block size becomes less than 32, reset dimensions to base values and select unshared code
   if (dB.x*dB.y < 32) {
      shared = 0;
      dB.x = basex;
      dB.y = basey;
   }

// Set the memory and grid dimensions
   tilelen=nn;
   mem = 0;
   if (shared > 0) mem = getsharedmemCP(dB,nn,mm);
   dG.x = (rkx+dB.x-1)/dB.x;  // x-dim over rank of 1st vector
   dG.y = (rky+dB.y-1)/dB.y;  // y-dim over rank of 2nd vector
   dG.z = (rkz+dB.z-1)/dB.z;  // z-dim over degrees-of-freedom
}

////////////////////////////////////////////////////////////////////////

__host__ void PVVkernel :: show()

////////////////////////////////////////////////////////////////////////
// Prints out parameters in PVV kernel

{
   cout << endl << " ---- PVV kernel ---- " << endl;
   cout << "n       = " << n << endl;
   cout << "m       = " << m << endl;
   cout << "rkx     = " << rkx << endl;
   cout << "rky     = " << rky << endl;
   cout << "tilelen = " << tilelen << endl;
   cout << "mem     = " << mem << endl;
   cout << "shared:   " << shared << endl;
   cout << "block dims : " << dB.x << " x " << dB.y << endl;
   cout << " grid dims : " << dG.x << " x " << dG.y << endl;
   cout << " -------------------- " << endl << endl;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__global__ void normbase_kernel(CPvec v)

////////////////////////////////////////////////////////////////////////
// Normalizes the 'base' of CP-vector, device kernel

{
   int j,k,n,kmod,s2;
   prec_typ norm1D,normND;

   int i  = blockIdx.x * blockDim.x + threadIdx.x;

   if (i<v.nrk){
      int ndof = v.ndof;
      int s1 = i*v.nrdim[ndof];

//    Loop over dimensions
      for (j=0,s2=s1,normND=1.0;j<ndof;j++){
          n=v.nbas[j];

//        1-D inner product
          for (k=0,norm1D=0.0;k<n;k++){
              kmod=s2+k;
              norm1D+=v.base[kmod]*v.base[kmod];
          }
          norm1D=sqrt(norm1D);

//        Normalization
          for (k=0;k<n;k++){
              kmod=s2+k;
              v.base[kmod]/=norm1D;
          }
          normND*=norm1D;
          s2+=n;
      }
      v.coef[i]*=normND;
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void normonebase_kernel(CPvec v, const int idof)

////////////////////////////////////////////////////////////////////////
// Normalizes the 'base' of CP-vector for mode 'idof', replacing the
// coef with the normalized value; device kernel

{
   int k;
   prec_typ norm1D,largest,tmp;

   int irk  = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk){
      int s = irk*v.nrdim[v.ndof] + v.nrdim[idof];
      int f = s + v.nbas[idof];

//    Get the largest entry
      largest=0.0;
      for (k=s;k<f;k++){
          largest=(largest > abs(v.base[k])) ? largest : abs(v.base[k]);
      }

//    1-D inner product
      for (k=s,norm1D=0.0;k<f;k++){
          tmp=v.base[k]/largest;
          norm1D+=tmp*tmp;
      }
      norm1D=largest*sqrt(abs(norm1D));

//    Normalize the base
      for (k=s;k<f;k++){
          v.base[k]/=norm1D;
      }

//    Replace the coefficient
      v.coef[irk]=norm1D;
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void normonebase_OLDkernel(CPvec v, const int idof)

////////////////////////////////////////////////////////////////////////
// Normalizes the 'base' of CP-vector for mode 'idof', replacing the
// coef with the normalized value; device kernel

{
   int k;
   prec_typ norm1D;

   int irk  = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk){
      int s = irk*v.nrdim[v.ndof] + v.nrdim[idof];
      int f = s + v.nbas[idof];

//    1-D inner product
      for (k=s,norm1D=0.0;k<f;k++){
          norm1D+=v.base[k]*v.base[k];
      }
      norm1D=sqrt(abs(norm1D));

//    Normalize the base
      for (k=s;k<f;k++){
          v.base[k]/=norm1D;
      }

//    Replace the coefficient
      v.coef[irk]=norm1D;
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void normalizecoef_kernel(CPvec v, const prec_typ *E)

////////////////////////////////////////////////////////////////////////
// Normalizes a CP-vec coef in device memory by dividing by sqrt(E),
// where 'E' is a normalization value computed elsewhere

{
// irk=rank index of v
   int irk = blockIdx.x * blockDim.x + threadIdx.x;

   if (irk < v.nrk) v.coef[irk] /= sqrt(abs(E[0]));
}

////////////////////////////////////////////////////////////////////////

__global__ void initps_kernel(CPvec v, CPvec w, Matrix PS, const bool multcoef)

////////////////////////////////////////////////////////////////////////
// Initializes 'PS' matrix of inner products to either ones or to
// products of coefficients, device kernel

{
   int vrk = v.nrk;
   int wrk = w.nrk;

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread computes a product and writes to PS
   if (i1 < vrk && i2 < wrk){
      prec_typ tmp=1.0;
      if (multcoef) tmp*=v.coef[i1]*w.coef[i2];
      SetMatrixelement(PS,i1,i2,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void multcoef_kernel(CPvec v, CPvec w, Matrix PS)

////////////////////////////////////////////////////////////////////////
// Multiplies 'PS' matrix of inner products by coefficients, device kernel

{
   int vrk = v.nrk;
   int wrk = w.nrk;

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread computes a product and writes to PS
   if (i1 < vrk && i2 < wrk){
      UpdateMatrixelement(PS,i1,i2,v.coef[i1]*w.coef[i2]);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constpsk_kernel(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset)

////////////////////////////////////////////////////////////////////////
// Updates PS matrix for the k-th dimension, device kernel
// Set reset=true to overwrite PS instead of update

{
   int j;
   int n=v.nbas[k];
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ tmp;

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Compute the inner product <v,w> if indices are in range; skip otherwise
   if (i1 < vrk && i2 < wrk){
      int i3 = v.nrdim[k];
      int nrdim = v.nrdim[v.ndof];
      int s1 = i1*nrdim+i3;
      int s2 = i2*nrdim+i3;
      for (j=0,tmp=0.0;j<n;j++){
          tmp+=v.base[s1+j]*w.base[s2+j];
      }
      if (reset){
         SetMatrixelement(PS,i1,i2,tmp);
      } else {
         UpdateMatrixelement(PS,i1,i2,tmp);
      }
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constpsk_shared_kernel(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset)

////////////////////////////////////////////////////////////////////////
// Updates 'PS' matrix of inner products for the k-th dimension, shared
// device kernel. Set reset=true to overwrite PS instead of update.

{
   int i3,j;
   int n=v.nbas[k];
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ tmp;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[n*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Load the base of vector 'v' into shared memory for the k-th dimension
   if (i1 < vrk){
      i3=threadIdx.y;
      while (i3 < n){
          tmp=GetCPelement(v,i1,k,i3);
          vs[i3+threadIdx.x*n]=tmp;
          i3+=blockDim.y;
      }
   }
   __syncthreads();

// Load vector 'w' also
   if (i2 < wrk){
      i3=threadIdx.x;
      while (i3 < n){
          tmp=GetCPelement(w,i2,k,i3);
          ws[i3+threadIdx.y*n]=tmp;
          i3+=blockDim.x;
      }
   }
   __syncthreads();

// Compute the inner product <v,w> if indices are in range
   if (i1 < vrk && i2 < wrk){
      int s1 = threadIdx.x*n;
      int s2 = threadIdx.y*n;
      for (j=0,tmp=0.0;j<n;j++){
          tmp+=vs[s1+j]*ws[s2+j];
      }
      if (reset){
         SetMatrixelement(PS,i1,i2,tmp);
      } else {
         UpdateMatrixelement(PS,i1,i2,tmp);
      }
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constpsk_shared_kernelA(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset)

////////////////////////////////////////////////////////////////////////
// Updates 'PS' matrix of inner products for the k-th dimension, shared
// device kernel. Set reset=true to overwrite PS instead of update.

{
   int i3,i4,j;
   int n=v.nbas[k];
   int ntile=(n+TILE_SIZE-1)/TILE_SIZE;
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ tmp,acc;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[TILE_SIZE*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;
   int s1 = threadIdx.x*TILE_SIZE;
   int s2 = threadIdx.y*TILE_SIZE;

   for (i4=0,acc=0.0;i4<ntile;i4++){

//    Load the base of vector 'v' into shared memory for the k-th dimension
      if (i1 < vrk){
         i3=threadIdx.y+i4*TILE_SIZE;
	 if (i3 < n){
             tmp=GetCPelement(v,i1,k,i3);
             vs[s1+threadIdx.y]=tmp;
	 }
      }
      __syncthreads();

//    Load vector 'w' also
      if (i2 < wrk){
	 i3=threadIdx.x+i4*TILE_SIZE;
	 if (i3 < n){
             tmp=GetCPelement(w,i2,k,i3);
             ws[s2+threadIdx.x]=tmp;
	 }
      }
      __syncthreads();

//    Compute the inner product <v,w>
      for (j=0;j<TILE_SIZE;j++){
	  i3=j+i4*TILE_SIZE;
	  if (i3 < n) acc+=vs[s1+j]*ws[s2+j];
      }
      __syncthreads();

   }

// Write the result to global array if indices are in range
   if (i1 < vrk && i2 < wrk){
      if (reset){
         SetMatrixelement(PS,i1,i2,acc);
      } else {
         UpdateMatrixelement(PS,i1,i2,acc);
      }
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constpsk_shared_kernelB(const CPvec v, const CPvec w, Matrix PS, const int k, const int tilelen, const bool reset)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, mode k, shared device kernel with tiling.

{
   int tos,mtos;
   int i,j,i3,s1,s2;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   int n=v.nbas[k];
   int ntile=(n+tilelen-1)/tilelen;
   int mos=v.nrdim[k]; // mode offset

   prec_typ acc;
   prec_typ tmp;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[tilelen*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Single-mode accumulator
   acc=0.0;

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){

       tos=itile*tilelen; // tile offset
       mtos=mos+tos;      // mode-&-tile offset

//     Load tile of v into shared memory
       for (i3=threadIdx.y;i3<tilelen;i3+=blockDim.y){
           tmp=0.0;
           j=tos+i3;
           if (i1 < vrk && j < n) tmp=v.base[i1*nrdim+mtos+i3];
           vs[i3+threadIdx.x*tilelen]=tmp;
       }
       __syncthreads();

//     Load tile of w into shared memory
       for (i3=threadIdx.x;i3<tilelen;i3+=blockDim.x){
           tmp=0.0;
           j=tos+i3;
           if (i2 < wrk && j < n) tmp=w.base[i2*nrdim+mtos+i3];
           ws[i3+threadIdx.y*tilelen]=tmp;
       }
       __syncthreads();

//     Compute the inner product: loop over tile length
       for (i=0,s1=threadIdx.x*tilelen,s2=threadIdx.y*tilelen;i<tilelen;i++,s1++,s2++){
           acc+=vs[s1]*ws[s2];
       }
       __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   if (i1 < vrk && i2 < wrk){
      if (reset){
         SetMatrixelement(PS,i1,i2,acc);
      } else {
         UpdateMatrixelement(PS,i1,i2,acc);
      }
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constpsk_registerized_kernel(const CPvec v, const CPvec w, Matrix PS, const int idof, const int tilelen, const bool reset)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different terms in
// the CP expansion, mode idof, shared device kernel with registerized loads.

{
   volatile int ilg,irg,icg,irow,icol,irk,iin;
   volatile int tos,mtos;
   int i,j,k;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   int n=v.nbas[idof];
   int ntile=(n+tilelen-1)/tilelen;
   int mos=v.nrdim[idof]; // mode offset
   prec_typ acc[WPT_X*WPT_Y];
   prec_typ wreg[WPT_Y];
   prec_typ vreg;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[tilelen*blockDim.x*WPT_X];

// Block offsets ({x,y} for rank-of-{v,w}, respectively)
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
       mtos=mos+tos;      // mode-&-tile offset

//     Load tile of v into shared memory
       for (ilg=0;ilg<nlgx;ilg++){
           irg=ilg/clgx; // row-group index (rank)
           icg=ilg%clgx; // col-group index (basis)
           irow=threadIdx.x+(irg+irsx*rlgx)*blockDim.x;
           icol=(threadIdx.y%tilelen)+icg*blockDim.y;
           irk=irow+bosx;
           iin=tos+icol;
           vreg=0.0;
           if (irk < vrk && iin < n) vreg=v.base[irk*nrdim+mtos+icol];
           vs[icol+irow*tilelen]=vreg;
       }
       __syncthreads();

//    Load tile of w into shared memory
      for (ilg=0;ilg<nlgy;ilg++){
          irg=ilg/clgy; // row-group index (rank)
          icg=ilg%clgy; // col-group index (basis)
          irow=threadIdx.y+(irg+irsy*rlgy)*blockDim.y;
          icol=(threadIdx.x%tilelen)+icg*blockDim.x;
          irk=irow+bosy;
          iin=tos+icol;
          vreg=0.0;
          if (irk < wrk && iin < n) vreg=w.base[irk*nrdim+mtos+icol];
          ws[icol+irow*tilelen]=vreg;
      }
      __syncthreads();

//    Compute the inner product: loop over tile length
      for (i=0;i<tilelen;i++){

//        Cache entries of w into registers
          for (k=0;k<WPT_Y;k++){
              icol=threadIdx.y+k*blockDim.y;
              wreg[k]=ws[i+icol*tilelen];
          }

//        Compute the inner product
          for (j=0;j<WPT_X;j++){
              irow=threadIdx.x+j*blockDim.x;
              vreg=vs[i+irow*tilelen]; // Cached v entry
              for (k=0;k<WPT_Y;k++){
                  acc[j+WPT_X*k]+=vreg*wreg[k];
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
	   if (irow < vrk && icol < wrk){
              if (reset){
	         SetMatrixelement(PS,irow,icol,acc[i]);
              } else {
                 UpdateMatrixelement(PS,irow,icol,acc[i]);
              }
	   }
       }
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constpstot_kernel(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, device kernel. Mode 'idof' is skipped, in
// case one desires the PS(!=idof) matrix; otherwise set idof=-1

{
   int j,k,n;
   prec_typ norm1D,normND;
   int vrk=v.nrk;
   int wrk=w.nrk;

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Compute the inner product <v,w> if indices are in range
   if (i1 < vrk && i2 < wrk){
      int s1 = i1*v.nrdim[v.ndof];
      int s2 = i2*w.nrdim[w.ndof];

      normND=1.0;
      if (multvcoef) normND*=v.coef[i1];
      if (multwcoef) normND*=w.coef[i2];

//    Loop over dimensions
      for (j=0;j<v.ndof;j++){

          n=v.nbas[j];

          if (j != idof) { // Skip mode 'idof'
//            1-D inner product
              for (k=0,norm1D=0.0;k<n;k++){
                  norm1D+=v.base[s1+k]*w.base[s2+k];
              }
              normND*=norm1D;
          }
          s1+=n;
          s2+=n;
      }
      SetMatrixelement(PS,i1,i2,normND);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constpstot_shared_kernel(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, shared device kernel. Mode 'idof' is 
// skipped, in case one desires the PS(!=idof) matrix; else set idof=-1

{
   int i3,j,k,n;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ norm1D,normND;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[nrdim*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Load the base of vector 'v' into shared memory
   if (i1 < vrk){
      i3=threadIdx.y;
      while (i3 < nrdim){
          vs[i3+threadIdx.x*nrdim]=v.base[i3+i1*nrdim];
          i3+=blockDim.y;
      }
   }
   __syncthreads();

// Load vector 'w' also
   if (i2 < wrk){
      i3=threadIdx.x;
      while (i3 < nrdim){
          ws[i3+threadIdx.y*nrdim]=w.base[i3+i2*nrdim];
          i3+=blockDim.x;
      }
   }
   __syncthreads();

// Compute the inner product <v,w> if indices are in range
   if (i1 < vrk && i2 < wrk){
      int s1 = threadIdx.x*nrdim;
      int s2 = threadIdx.y*nrdim;

      normND=1.0;
      if (multvcoef) normND*=v.coef[i1];
      if (multwcoef) normND*=w.coef[i2];

//    Loop over dimensions
      for (j=0;j<ndof;j++){
          n=v.nbas[j];
          if (j != idof) { // Skip mode 'idof'
//            1-D inner product
              for (k=0,norm1D=0.0;k<n;k++){
                  norm1D+=vs[s1+k]*ws[s2+k];
              }
              normND*=norm1D;
          }
          s1+=n;
          s2+=n;
      }
      SetMatrixelement(PS,i1,i2,normND);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void constpstot_shared_kernelA(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, shared device kernel. Mode 'idof' is 
// skipped, in case one desires the PS(!=idof) matrix; else set idof=-1
// This version builds PS mode-by mode sequentially
{
   int i3,j,k,n;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   int nmax=v.nmax;
   prec_typ norm1D,normND;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[nmax*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

   normND=1.0;
   if (i1 < vrk && multvcoef) normND*=v.coef[i1];
   if (i2 < wrk && multwcoef) normND*=w.coef[i2];

// Loop over dimensions
   for (j=0;j<ndof;j++){
       n=v.nbas[j];
       int s=v.nrdim[j];

//     Load the base of j-th mode of vector 'v' into shared memory
       if (i1 < vrk){
          i3=threadIdx.y;
          while (i3 < n){
              vs[i3+threadIdx.x*nmax]=v.base[i3+s+i1*nrdim];
              i3+=blockDim.y;
          }
       }
       __syncthreads();

//     Load vector 'w' also
       if (i2 < wrk){
          i3=threadIdx.x;
          while (i3 < n){
              ws[i3+threadIdx.y*nmax]=w.base[i3+s+i2*nrdim];
              i3+=blockDim.x;
          }
       }
       __syncthreads();

//     Compute j-th mode inner product <v,w> if indices are in range
       if (i1 < vrk && i2 < wrk){
          int s1 = threadIdx.x*nmax;
          int s2 = threadIdx.y*nmax;

          if (j != idof) { // Skip mode 'idof'
//           1-D inner product
             for (k=0,norm1D=0.0;k<n;k++){
                 norm1D+=vs[s1+k]*ws[s2+k];
             }
             normND*=norm1D;
          }
       }
   }

// Write final result to global memory
   if (i1 < vrk && i2 < wrk) SetMatrixelement(PS,i1,i2,normND);

}

////////////////////////////////////////////////////////////////////////

__global__ void constpstot_shared_kernelB(const CPvec v, const CPvec w, Matrix PS, const int idof, const int tilelen, const bool multvcoef, const bool multwcoef)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, shared device kernel with tilint.
// Mode 'idof' is skipped, in case one desires the PS(!=idof) matrix; else set idof=-1

{
   int n,ntile,mos,tos,mtos;
   int i,j,i3,s1,s2;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ acc1;
   prec_typ accD;
   prec_typ tmp;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[tilelen*blockDim.x];

// i1=row index, i2=col index, of PS
   int i1 = blockIdx.x * blockDim.x + threadIdx.x;
   int i2 = blockIdx.y * blockDim.y + threadIdx.y;

// Multi-mode accumulator
   accD=1.0;

// Loop over modes, skipping idof
   for (int d=0;d<ndof;d++){

       if (d != idof){
          n=v.nbas[d];
          ntile=(n+tilelen-1)/tilelen;
          mos=v.nrdim[d]; // mode offset

//        Single-mode accumulator
          acc1=0.0;

//        Loop over tiles
          for (int itile=0;itile<ntile;itile++){

              tos=itile*tilelen; // tile offset
              mtos=mos+tos;      // mode-&-tile offset

//            Load tile of v into shared memory
              for (i3=threadIdx.y;i3<tilelen;i3+=blockDim.y){
                  tmp=0.0;
                  j=tos+i3;
                  if (i1 < vrk && j < n) tmp=v.base[i1*nrdim+mtos+i3];
                  vs[i3+threadIdx.x*tilelen]=tmp;
              }
              __syncthreads();

//           Load tile of w into shared memory
              for (i3=threadIdx.x;i3<tilelen;i3+=blockDim.x){
                  tmp=0.0;
                  j=tos+i3;
                  if (i2 < wrk && j < n) tmp=w.base[i2*nrdim+mtos+i3];
                  ws[i3+threadIdx.y*tilelen]=tmp;
              }
             __syncthreads();

//           Compute the inner product: loop over tile length
             for (i=0,s1=threadIdx.x*tilelen,s2=threadIdx.y*tilelen;i<tilelen;i++,s1++,s2++){
		 acc1+=vs[s1]*ws[s2];
             }
             __syncthreads();

          } // Loop over tiles (itile)

//        Multiply the 1-mode product with the multi-mode accumulated product
          accD*=acc1;

       } // if not skipped mode
   } // Loop over DOF (d)

// Write the final result to global memory
   if (i1 < vrk && i2 < wrk){
      if (multvcoef) accD*=v.coef[i1];
      if (multwcoef) accD*=w.coef[i2];
      SetMatrixelement(PS,i1,i2,accD);
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void constpstot_registerized_kernel(const CPvec v, const CPvec w, Matrix PS, const int idof, const int tilelen, const bool multvcoef, const bool multwcoef)

////////////////////////////////////////////////////////////////////////
// Constructs 'PS' matrix, the matrix of inner products of different
// terms in the CP expansion, shared device kernel with registerized loads.
// Mode 'idof' is skipped, in case one desires the PS(!=idof) matrix; else set idof=-1

{
   volatile int ilg,irg,icg,irow,icol,irk,iin;
   volatile int n,ntile,mos,tos,mtos;
   int i,j,k;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   int wrk=w.nrk;
   prec_typ acc1[WPT_X*WPT_Y];
   prec_typ accD[WPT_X*WPT_Y];
   prec_typ wreg[WPT_Y];
   prec_typ vreg;

   extern __shared__ prec_typ vws[];
   prec_typ* vs = (prec_typ*) vws;
   prec_typ* ws = (prec_typ*) &vs[tilelen*blockDim.x*WPT_X];

// Block offsets ({x,y} for rank-of-{v,w}, respectively)
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

// Multi-mode product registers
   for (i=0;i<WPT_X*WPT_Y;i++) accD[i]=1.0;

// Loop over modes, skipping idof
   for (int d=0;d<ndof;d++){

       if (d != idof){
          n=v.nbas[d];
          ntile=(n+tilelen-1)/tilelen;
          mos=v.nrdim[d]; // mode offset

//        Single-mode product registers
          for (i=0;i<WPT_X*WPT_Y;i++) acc1[i]=0.0;

//        Loop over tiles
          for (int itile=0;itile<ntile;itile++){

              tos=itile*tilelen; // tile offset
              mtos=mos+tos;      // mode-&-tile offset

//            Load tile of v into shared memory
              for (ilg=0;ilg<nlgx;ilg++){
                  irg=ilg/clgx; // row-group index (rank)
                  icg=ilg%clgx; // col-group index (basis)
                  irow=threadIdx.x+(irg+irsx*rlgx)*blockDim.x;
                  icol=(threadIdx.y%tilelen)+icg*blockDim.y;
                  irk=irow+bosx;
                  iin=tos+icol;
                  vreg=0.0;
                  if (irk < vrk && iin < n) vreg=v.base[irk*nrdim+mtos+icol];
                  vs[icol+irow*tilelen]=vreg;
              }
              __syncthreads();

//           Load tile of w into shared memory
             for (ilg=0;ilg<nlgy;ilg++){
                 irg=ilg/clgy; // row-group index (rank)
                 icg=ilg%clgy; // col-group index (basis)
                 irow=threadIdx.y+(irg+irsy*rlgy)*blockDim.y;
                 icol=(threadIdx.x%tilelen)+icg*blockDim.x;
                 irk=irow+bosy;
                 iin=tos+icol;
                 vreg=0.0;
                 if (irk < wrk && iin < n) vreg=w.base[irk*nrdim+mtos+icol];
                 ws[icol+irow*tilelen]=vreg;
             }
             __syncthreads();

//           Compute the inner product: loop over tile length
             for (i=0;i<tilelen;i++){

//               Cache entries of w into registers
                 for (k=0;k<WPT_Y;k++){
                     icol=threadIdx.y+k*blockDim.y;
                     wreg[k]=ws[i+icol*tilelen];
                 }

//               Compute the inner product
	         for (j=0;j<WPT_X;j++){
                     irow=threadIdx.x+j*blockDim.x;
		     vreg=vs[i+irow*tilelen]; // Cached v entry
                     for (k=0;k<WPT_Y;k++){
                         acc1[j+WPT_X*k]+=vreg*wreg[k];
                     }
                 }

             }
             __syncthreads();

          } // Loop over tiles (itile)

//        Multiply the 1-mode product with the multi-mode accumulated product
          for (i=0;i<WPT_X*WPT_Y;i++) accD[i]*=acc1[i];

       } // if not skipped mode

   } // Loop over DOF (d)

// Cache coefs of w into registers
   if (multwcoef) {
      for (k=0;k<WPT_Y;k++){
          icol=bosy+threadIdx.y+k*blockDim.y;
          if (icol < wrk) wreg[k]=w.coef[icol];
      }
   }

// Write the final result to global memory
   for (j=0;j<WPT_X;j++){
       irow=bosx+threadIdx.x+j*blockDim.x;
       if (irow < vrk) vreg=v.coef[irow];
       for (k=0;k<WPT_Y;k++){
           icol=bosy+threadIdx.y+k*blockDim.y;
	   i=j+WPT_X*k;
	   if (irow < vrk && icol < wrk){
              if (multvcoef) accD[i]*=vreg;
              if (multwcoef) accD[i]*=wreg[k];
              SetMatrixelement(PS,irow,icol,accD[i]);
	   }
       }
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void sumoverrow_kernel(prec_typ *v, const Matrix M)

////////////////////////////////////////////////////////////////////////
// Selects a matrix row and sums the column elements, device kernel

{
   int i;
   int nr = M.nr;
   int row = blockIdx.x * blockDim.x + threadIdx.x;

   if (row < nr){
      for (i=0,v[row]=0.0;i<M.nc;i++){
          v[row]+=M.mat[row+i*nr];
      }
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void sumovercol_kernel(prec_typ *v, const Matrix M)

////////////////////////////////////////////////////////////////////////
// Selects a matrix column and sums the row elements, device kernel

{
   int i;
   int nr = M.nr;
   int col = blockIdx.x * blockDim.x + threadIdx.x;
   int ist = col*nr;

   if (col < M.nc){
      for (i=0,v[col]=0.0;i<nr;i++){
          v[col]+=M.mat[ist+i];
      }
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void sumarray_kernel(prec_typ *sum, const prec_typ *v, const int n)

////////////////////////////////////////////////////////////////////////
// Sums elements in array v, device kernel

{
   int i;

   for (i=0,sum[0]=0.0;i<n;i++){
       sum[0]+=v[i];
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-specific functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void constpsk_wrapper(const PVVkernel &PV, const CPvec &v, const CPvec &w, Matrix &PS, const int &idof, const bool &reset)

////////////////////////////////////////////////////////////////////////
// Wrapper for shared/unshared calls to constpsk kernel

{
   if (PV.shared == 2) {
      constpsk_registerized_kernel<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,PV.tilelen,reset);
   } else if (PV.shared == 1) {
//      constpsk_shared_kernel<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,reset);
//      constpsk_shared_kernelA<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,reset);
      constpsk_shared_kernelB<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,PV.tilelen,reset);
   } else {
      constpsk_kernel<<<PV.dG,PV.dB>>>(v,w,PS,idof,reset);
   }
}

////////////////////////////////////////////////////////////////////////

__host__ void constpstot_wrapper(const PVVkernel &PV, const CPvec &v, const CPvec &w, Matrix &PS, const int &idof, const bool &multvcoef, const bool &multwcoef)

////////////////////////////////////////////////////////////////////////
// Wrapper for shared/unshared calls to constpstot kernel

{
   if (PV.shared == 2) {
       constpstot_registerized_kernel<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,PV.tilelen,multvcoef,multwcoef);
   } else if (PV.shared == 1) {
//      constpstot_shared_kernel<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,multvcoef,multwcoef);
//      constpstot_shared_kernelA<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,multvcoef,multwcoef);
      constpstot_shared_kernelB<<<PV.dG,PV.dB,PV.mem>>>(v,w,PS,idof,PV.tilelen,multvcoef,multwcoef);
   } else {
      constpstot_kernel<<<PV.dG,PV.dB>>>(v,w,PS,idof,multvcoef,multwcoef);
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
