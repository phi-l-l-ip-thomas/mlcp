//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CUDA-driven "intertwining" solver, MVP portion
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains the CUDA version of the intertwining code, matrix-vector product part

#include <iostream>
#include <stdio.h>
#include <omp.h>
#include "global.h"
#include "utils.h"
#include "matrix.h"
#include "cpvector.h"
#include "hoplist.h"
#include "modvecvec.h"
#include "linalg.h"
#include "als.h"
#include "intertw_mvp.h"

using namespace std;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__global__ void mvp_kernel(const HopList HL, const int iterm, const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Applies CP-format matrix-vector product, Hv=w, device kernel

{
   int i,ir,ic;
   prec_typ tmp;

// irk=rank index, ibas=row-of-H basis index, idof=dimension index
   int irk  = blockIdx.x * blockDim.x + threadIdx.x;
   int ibas = blockIdx.y * blockDim.y + threadIdx.y;
   int idof = blockIdx.z * blockDim.z + threadIdx.z;

// CP sizes
   int n = v.nbas[idof];
   int rk = v.nrk;

// Hamiltonian components
// hsft tells where in HL to find the appropriate h matrix
   int hsft = HL.mptr[HL.opid[iterm*v.ndof+idof]-1]-1;
   prec_typ alpha = 1.0;
   if (idof == 0) alpha = HL.coef[iterm];

// Compute an element of w if indices are in range
   if (irk < rk && ibas < n){
      for (i=0,tmp=0.0;i<n;i++){
          ir=min(i,ibas);  // ir,ic are row and column of the h-matrix
          ic=max(i,ibas);
          tmp+=v.base[irk*v.nrdim[v.ndof]+v.nrdim[idof]+i]*HL.mat[hsft+ir+ic*(ic+1)/2];

      }
      tmp*=alpha;
      SetCPelement(w,irk,idof,ibas,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void mvp_shared_kernel(const HopList HL, const int iterm, const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Applies CP-format matrix-vector product, Hv=w, shared device kernel
 
{
   int i,ir,ic;
   prec_typ tmp;

// irk=rank index, ibas=row-of-H basis index, idof=dimension index
   int irk  = blockIdx.x * blockDim.x + threadIdx.x;
   int ibas = blockIdx.y * blockDim.y + threadIdx.y;
   int idof = blockIdx.z * blockDim.z + threadIdx.z;

// CP sizes
   int n = v.nbas[idof];
   int rk = v.nrk;

// Hamiltonian components
// hsft tells where in HL to find the appropriate h matrix
   int hsft = HL.mptr[HL.opid[iterm*v.ndof+idof]-1]-1;
   prec_typ alpha = 1.0;
   if (idof == 0) alpha = HL.coef[iterm];

// Allocatable arrays
   extern __shared__ prec_typ vHs[];
   prec_typ* vs = (prec_typ*) vHs;
   prec_typ* Hs = (prec_typ*) &vs[n*blockDim.x];

// Load the base of the idof dimension of vector 'v' into shared memory
   if (irk < rk){
      i=threadIdx.y;
      while (i < n){
          tmp=GetCPelement(v,irk,idof,i);
          vs[i+threadIdx.x*n]=tmp;
          i+=blockDim.y;
      }
   }
   __syncthreads();

// Load rows of H into shared memory
   if (ibas < n){
      i=threadIdx.x;
      while (i < n){
          ir=min(i,ibas);  // ir,ic are row and column of the h-matrix
          ic=max(i,ibas);
          Hs[i+threadIdx.y*n]=HL.mat[hsft+ir+ic*(ic+1)/2];
          i+=blockDim.x;
      }
   }
   __syncthreads();

// Compute an element of w if indices are in range
   if (irk < rk && ibas < n){
      int s1 = threadIdx.x*n;
      int s2 = threadIdx.y*n;
      for (i=0,tmp=0.0;i<n;i++){
          tmp+=vs[s1+i]*Hs[s2+i];
      }
      tmp*=alpha;
      SetCPelement(w,irk,idof,ibas,tmp);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void mvp_shared_kernelB(const HopList HL, const int iterm, const int tilelen, const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Applies CP-format matrix-vector product, Hv=w, shared device kernel
 
{
   int i,j,ir,ic,tos,mtos;
   prec_typ tmp,acc;

// irk=rank index, ibas=row-of-H basis index, idof=dimension index
   int irk  = blockIdx.x * blockDim.x + threadIdx.x;
   int ibas = blockIdx.y * blockDim.y + threadIdx.y;
   int idof = blockIdx.z * blockDim.z + threadIdx.z;

// CP sizes
   int n = v.nbas[idof];
   int mos=v.nrdim[idof]; // mode offset
   int ntile=(n+tilelen-1)/tilelen;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int rk = v.nrk;

// Hamiltonian components
// hsft tells where in HL to find the appropriate h matrix
   int hsft = HL.mptr[HL.opid[iterm*v.ndof+idof]-1]-1;
   prec_typ alpha = 1.0;
   if (idof == 0) alpha = HL.coef[iterm];

// Allocatable arrays
   extern __shared__ prec_typ vHs[];
   prec_typ* vs = (prec_typ*) vHs;
   prec_typ* Hs = (prec_typ*) &vs[tilelen*blockDim.x];

// Single-mode accumulator
   acc=0.0;

// Loop over tiles
   for (int itile=0;itile<ntile;itile++){  

       tos=itile*tilelen; // tile offset
       mtos=mos+tos;      // mode-&-tile offset

//     Load tile of v into shared memory
       for (i=threadIdx.y;i<tilelen;i+=blockDim.y){
           tmp=0.0;
	   j=tos+i;
           if (irk < rk && j < n) tmp=v.base[irk*nrdim+mtos+i];
           vs[i+threadIdx.x*tilelen]=tmp;
       }
       __syncthreads();

//     Load rows of H into shared memory
       for (i=threadIdx.x;i<tilelen;i+=blockDim.x){
           tmp=0.0;
	   j=tos+i;
           if (ibas < n && j < n) {
              ir=min(j,ibas);  // ir,ic are row and column of the h-matrix
              ic=max(j,ibas);
              tmp=HL.mat[hsft+ir+ic*(ic+1)/2];
	   }
	   Hs[i+threadIdx.y*tilelen]=tmp;
       }
       __syncthreads();

//     Compute the inner product: loop over tile length
       for (i=0,ir=threadIdx.x*tilelen,ic=threadIdx.y*tilelen;i<tilelen;i++,ir++,ic++){
           acc+=vs[ir]*Hs[ic];
       }
       __syncthreads();

   }

// Copy element of w to global if indices are in range
   if (irk < rk && ibas < n){
      acc*=alpha;
      SetCPelement(w,irk,idof,ibas,acc);
   }

}

////////////////////////////////////////////////////////////////////////

__global__ void mvp_registerized_kernel(const HopList HL, const int iterm, const int tilelen, const CPvec v, CPvec w)

////////////////////////////////////////////////////////////////////////
// Applies CP-format matrix-vector product, Hv=w, shared device kernel with registerized loads.

{
   volatile int ilg,irg,icg,irow,icol,irk,iin;
   volatile int tos,mtos;
   int i,j,k,hr,hc;
   int ndof=v.ndof;
   int nrdim=v.nrdim[ndof];
   int vrk=v.nrk;
   prec_typ acc[WPT_X*WPT_Y];
   prec_typ Hreg[WPT_Y];
   prec_typ vreg;

   extern __shared__ prec_typ vHs[];
   prec_typ* vs = (prec_typ*) vHs;
   prec_typ* Hs = (prec_typ*) &vs[tilelen*blockDim.x*WPT_X];

// Hamiltonian components
// hsft tells where in HL to find the appropriate h matrix
   int idof = blockIdx.z * blockDim.z + threadIdx.z;
   int hsft = HL.mptr[HL.opid[iterm*ndof+idof]-1]-1;
   prec_typ alpha = 1.0;
   if (idof == 0) alpha = HL.coef[iterm];

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

// Mode-specific values
   int n=v.nbas[idof];
   int ntile=(n+tilelen-1)/tilelen;
   int mos=v.nrdim[idof]; // mode offset

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

//    Load tile of H into shared memory
      for (ilg=0;ilg<nlgy;ilg++){
          irg=ilg/clgy; // row-group index (basis row)
          icg=ilg%clgy; // col-group index (basis col)
          irow=threadIdx.y+(irg+irsy*rlgy)*blockDim.y;
          icol=(threadIdx.x%tilelen)+icg*blockDim.x;
          irk=irow+bosy; // actually a row-of-H, not a rank
          iin=tos+icol;
          vreg=0.0;
          if (irk < n && iin < n){
             hr=min(irk,iin);
             hc=max(irk,iin);
	     vreg=HL.mat[hsft+hr+hc*(hc+1)/2];
          }
          Hs[icol+irow*tilelen]=vreg;
      }
      __syncthreads();

//    Compute the inner product: loop over tile length
      for (i=0;i<tilelen;i++){

//        Cache entries of H into registers
          for (k=0;k<WPT_Y;k++){
              icol=threadIdx.y+k*blockDim.y;
              Hreg[k]=Hs[i+icol*tilelen];
          }

//        Compute the inner product
          for (j=0;j<WPT_X;j++){
              irow=threadIdx.x+j*blockDim.x;
              vreg=vs[i+irow*tilelen]; // Cached v entry
              for (k=0;k<WPT_Y;k++){
                  acc[j+WPT_X*k]+=vreg*Hreg[k];
              }
          }
      }
      __syncthreads();

   } // Loop over tiles (itile)

// Write the final result to global memory
   for (j=0;j<WPT_X;j++){
       irow=bosx+threadIdx.x+j*blockDim.x;
       if (irow < vrk) vreg=v.coef[irow];
       for (k=0;k<WPT_Y;k++){
           icol=bosy+threadIdx.y+k*blockDim.y;
	   i=j+WPT_X*k;
	   if (irow < vrk && icol < n){
              acc[i]*=alpha;
	      SetCPelement(w,irow,idof,icol,acc[i]);
	   }
       }
   }

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void mvp_wrapper(const PVVkernel &PV, const HopList &HL, const int &iterm, const CPvec &v, CPvec &w)

////////////////////////////////////////////////////////////////////////
// Wrapper for shared/unshared calls to constbjk kernel

{
   if (PV.shared == 2) {
      mvp_registerized_kernel<<<PV.dG,PV.dB,PV.mem>>>(HL,iterm,PV.tilelen,v,w);
   } else if (PV.shared == 1) {
//      mvp_shared_kernel<<<PV.dG,PV.dB,PV.mem>>>(HL,iterm,v,w);
      mvp_shared_kernelB<<<PV.dG,PV.dB,PV.mem>>>(HL,iterm,PV.tilelen,v,w);
   } else {
      mvp_kernel<<<PV.dG,PV.dB>>>(HL,iterm,v,w);
   }
}

////////////////////////////////////////////////////////////////////////

__host__ void mvp_build_linear_system(const PVVkernel mvp, const PVVkernel pst, const PVVkernel bjk, const HopList &HL_d, const CPvec &Fo_d, CPvec &F_d, CPvec &G_d, Matrix &BB_d, Matrix &bjk_d, const int &idof, const int &ni, const int &ishift, const prec_typ &Eshift, const bool recalcB)

////////////////////////////////////////////////////////////////////////
// Builds the BB, bjk matrices required for solving the linear systems
// for the intertwined power method.

{
   int k;

// Initialize values
   int T = HL_d.nrk;
   int rk = F_d.nrk;
   int nD = F_d.nD;
   prec_typ beta = 0.0;

// Initialize CUDA types
   dim3 dB2(BLOCK_SIZE_2D,BLOCK_SIZE_2D);
   dim3 dB1(BLOCK_SIZE_1D);
   dim3 dG1((rk+dB1.x-1)/dB1.x);                             // rk x 1
   dim3 dimGridbjk((rk+dB2.x-1)/dB2.x,(ni+dB2.y-1)/dB2.y);   // rk x n
   dim3 dimGridCopy((rk+dB2.x-1)/dB2.x,(nD+dB2.y-1)/dB2.y);  // rk x nD

// Error checking
   if (Fo_d.nrk != rk){
      cout << "rk, rkold = " << rk << ", " << Fo_d.nrk << " must be equal" << endl;
      throwERROR("error in mvp_build_linear_system()");
   }

// Initialize bjk on device
   fillmatrix_kernel<<<dimGridbjk,dB2>>>(bjk_d,beta);

// Copy coefficients of F to G on device
   copycpcoef_kernel<<<dG1,dB1>>>(Fo_d,G_d);

   for (k=0;k<T;k++){ // Loop over terms in H // TEST
//     Apply matrix-vector product H*F = G
       mvp_wrapper(mvp,HL_d,k,Fo_d,G_d);

//     Construct PS(!=j) matrix
       constpstot_wrapper(pst,G_d,F_d,BB_d,idof,false,false);

//     Add terms to bjk matrix
       constbjk_wrapper(bjk,G_d,BB_d,bjk_d,idof);
   }

// Copy the base of Fo_d -> G_d for constructing BB and applying the energy shift, if applicable
   copycpbase_kernel<<<dimGridCopy,dB2>>>(Fo_d,G_d);

// Construct BB(!=j) matrix, which here is also the PS(!=j) matrix for the energy shift
   constpstot_wrapper(pst,G_d,F_d,BB_d,idof,false,false);

// Add terms to bjk matrix
   if (ishift != 0){
      scalecpbase_kernel<<<dimGridbjk,dB2>>>(G_d,idof,-Eshift);   // Apply energy shift basis of mode 'idof' in G
      constbjk_wrapper(bjk,G_d,BB_d,bjk_d,idof);
      if (ishift < 0){  // (E1-H)*F instead of (H-E1)*F: change sign of bjk
         scalematrix_kernel<<<dimGridbjk,dB2>>>(bjk_d,-1.0);
      }
   }

// In the case Fo_d != F_d, BB_d must be recalculated since it is <F_d,F_d>
   if (recalcB){
//    Copy the base of F_d -> G_d
      copycpbase_kernel<<<dimGridCopy,dB2>>>(F_d,G_d);

//    Construct BB(!=j) matrix, this time different from PS(!=j)
      constpstot_wrapper(pst,G_d,F_d,BB_d,idof,false,false);
   }
}

////////////////////////////////////////////////////////////////////////

__host__ void POW_ALS(const HopList &HL_d, intwvars &IV, const int &npow, const int &ishift, const prec_typ &Eshift)

////////////////////////////////////////////////////////////////////////
// Intertwined power method. Given CP-format Hamiltonian H and vector F
// F is iteratively refined using the power method (H - Eshift*I)*F.

{
   if (npow == 0) return;

   int i,j,n;
   int rk = IV.rk;
   int ndof = IV.F_d.ndof;

// Compute the regularization penalty from the coefs of F
   setpenalty_kernel<<<IV.dG,IV.dB>>>(IV.pen_d,IV.F_d.coef,rk);

// Loop over intertwining interations
   for (i=0;i<npow;i++){

//     Loop over DOF
       for (j=0;j<ndof;j++){
           n = IV.nbas[j];

//         Calculate BB and b_j_k
           mvp_build_linear_system(IV.mvp,IV.pst,IV.bjk,HL_d,IV.F_d,IV.F_d,IV.G_d,IV.BB_d,IV.bjk_d,j,n,ishift,Eshift,false);

//         Add penalty to avoid ill-conditioning
           addpenalty_kernel<<<IV.dG1,IV.dB1>>>(IV.BB_d,IV.pen_d);

//         Solve linear system B*c_j_k = b_j_k
           solve_linsys(IV.BB_d,IV.bjk_d,IV.WORK,IV.IPV,IV.INFO,rk,n);

//         Construct improved F and normalize base
           dim3 dimGridbjk((rk+IV.dB2.x-1)/IV.dB2.x,(n+IV.dB2.y-1)/IV.dB2.y); // rk x n
           replacevwithsolns_kernel<<<dimGridbjk,IV.dB2>>>(IV.F_d,IV.bjk_d,j);
           normonebase_kernel<<<IV.dG1,IV.dB1>>>(IV.F_d,j);
       }

//     Normalize coefficients at the end of the sweep
       constpstot_wrapper(IV.pst,IV.F_d,IV.F_d,IV.BB_d,-1,true,true); // Reconstruct BB_d
       sumoverrow_kernel<<<IV.dG1,IV.dB1>>>(IV.tmp_d,IV.BB_d);   // Sum over rows of BB_d
       sumarray_kernel<<<IV.dG,IV.dB>>>(IV.sum_d,IV.tmp_d,rk);   // Accumulate sum over column
       normalizecoef_kernel<<<IV.dG1,IV.dB1>>>(IV.F_d,IV.sum_d); // Normalize: divide by sqrt(sum_d)
   }

}

////////////////////////////////////////////////////////////////////////

__host__ void PRODHV_ALS(const HopList &HL_d, intwvars &IV, const int &npow, const int &ishift, const prec_typ &Eshift)

////////////////////////////////////////////////////////////////////////
// Intertwined matrix-vector product. Given CP-format Hamiltonian H and 
// vector F, F is computed/refined for a single matrix-vector product,
// (H - Eshift*I)*F.

{
   if (npow == 0) return;

   int i,j,n;
   int rk = IV.rk;
   int ndof = IV.F_d.ndof;
   int nD = IV.nD;

// Initialize CUDA types
   dim3 dimGridCopy((rk+IV.dB2.x-1)/IV.dB2.x,(nD+IV.dB2.y-1)/IV.dB2.y); // rk x nD

// Copy Fo -> F on device (both coefs and base)
   copycpcoef_kernel<<<IV.dG1,IV.dB1>>>(IV.Fo_d,IV.F_d);
   copycpbase_kernel<<<dimGridCopy,IV.dB2>>>(IV.Fo_d,IV.F_d);

// Compute the regularization penalty from the coefs of F
   setpenalty_kernel<<<IV.dG,IV.dB>>>(IV.pen_d,IV.Fo_d.coef,rk);

// Loop over intertwining interations
   for (i=0;i<npow;i++){

//     Loop over DOF
       for (j=0;j<ndof;j++){
           n = IV.nbas[j];

//         Calculate BB and b_j_k
           mvp_build_linear_system(IV.mvp,IV.pst,IV.bjk,HL_d,IV.Fo_d,IV.F_d,IV.G_d,IV.BB_d,IV.bjk_d,j,n,ishift,Eshift,true);

//         Add penalty to avoid ill-conditioning
           addpenalty_kernel<<<IV.dG1,IV.dB1>>>(IV.BB_d,IV.pen_d);

//         Solve linear system B*c_j_k = b_j_k
           solve_linsys(IV.BB_d,IV.bjk_d,IV.WORK,IV.IPV,IV.INFO,rk,n);

//         Construct improved F and normalize base
           dim3 dimGridbjk((rk+IV.dB2.x-1)/IV.dB2.x,(n+IV.dB2.y-1)/IV.dB2.y); // rk x n
           replacevwithsolns_kernel<<<dimGridbjk,IV.dB2>>>(IV.F_d,IV.bjk_d,j);
           normonebase_kernel<<<IV.dG1,IV.dB1>>>(IV.F_d,j);
       }
   }

}

////////////////////////////////////////////////////////////////////////

__host__ void MatrixVectorProducts(const HopList &HL_d, CPvec &Q_h, const int &ncpu, const int &npow, const int &nbloc, const int &rk, const int &ishift, const prec_typ &Eshift)

////////////////////////////////////////////////////////////////////////
// Performs the intertwined power method, matrix-vector product portion,
// on a block of vectors

{
   int i,j;

// Initialize MVP containers
   intwvars IV[ncpu];
   for (j=0;j<ncpu;j++){
       IV[j].newmvp(Q_h,rk);
   }

// Perform intertwined power iterations
#pragma omp parallel shared(IV) private(i,j)
{
#pragma omp for
   for (i=0;i<nbloc;i++){
       j=omp_get_thread_num();
//     Send vector to device, compute MVPs, copy back to host
       IV[j].send(Q_h,i);
       POW_ALS(HL_d,IV[j],npow,ishift,Eshift);
       IV[j].get(Q_h,i);
   }
}

// Dispose MVP containers
   for (j=0;j<ncpu;j++){
       IV[j].flush();
   }
}

////////////////////////////////////////////////////////////////////////

__host__ void GenEigvalProblem(const HopList &HL_d, const CPvec &Q_h, CPvec &Q_d, Matrix &QHQ_d, prec_typ *eigv_h, const int &ncpu, const int &nals, const int &nbloc, const int &rk)

////////////////////////////////////////////////////////////////////////
// Solves the generalized eigenvalue problem on device

{
   int i,j,k;
   prec_typ *eigv_d;
   Matrix S_d;
   CPvec G_d[ncpu];

// Arrays for the generalized eigenvalue problem
   NewMatrixonDevice(S_d,nbloc,nbloc);
   cudaMalloc(&eigv_d,nbloc*sizeof(prec_typ));

// Initialize Generalized Eigenvalue Problem containers
   intwvars IV[ncpu];
   for (k=0;k<ncpu;k++){
       IV[k].newgep(Q_h,rk,nbloc);
   }

// Loop over vectors in the block
#pragma omp parallel shared(IV,QHQ_d,S_d,G_d) private(i,j,k)
{
#pragma omp for
   for (i=0;i<nbloc;i++){
       k=omp_get_thread_num();

//     Select a vector in the block; designate as "Fo"
       IV[k].point(Q_d,i);

//     Compute H*F
       PRODHV_ALS(HL_d,IV[k],nals,0,0.0);

//     Inner products: <G,F> and <G,H*F>
       for (j=0;j<i+1;j++){
           PointCPvecRtermsonGPU(Q_d,G_d[k],j*rk,rk);

//         Build QHQ(i,j)=<Q(j),H*Q(i)>
           constpstot_wrapper(IV[k].pst,G_d[k],IV[k].F_d,IV[k].BB_d,-1,true,true);  // Build BB_d
           sumoverrow_kernel<<<IV[k].dG1,IV[k].dB1>>>(IV[k].tmp_d,IV[k].BB_d); // Sum over rows of BB_d
           sumarray_kernel<<<IV[k].dG,IV[k].dB>>>(&QHQ_d.mat[i+j*nbloc],IV[k].tmp_d,rk); // Accumulate sum in QHQ

//         Build S(i,j)=<Q(j),Q(i)>
           constpstot_wrapper(IV[k].pst,G_d[k],IV[k].Fo_d,IV[k].BB_d,-1,true,true); // Build BB_d
           sumoverrow_kernel<<<IV[k].dG1,IV[k].dB1>>>(IV[k].tmp_d,IV[k].BB_d); // Sum over rows of BB_d
           sumarray_kernel<<<IV[k].dG,IV[k].dB>>>(&S_d.mat[i+j*nbloc],IV[k].tmp_d,rk); // Accumulate sum in S
       }
   }
}

// Dispose Generalized Eigenvalue Problem containers
   for (k=0;k<ncpu;k++){
       IV[k].flush();
   }

// Solve the generalized eigenvalue problem
   solve_geneigval(QHQ_d,S_d,eigv_d,nbloc);

// Copy eigenvalues host <- device, overwriting eigv_h on host
   cudaMemcpy(eigv_h,eigv_d,nbloc*sizeof(prec_typ),cudaMemcpyDeviceToHost);

// Cleanup
   cudaFree(eigv_d);
   FlushMatrixonDevice(S_d);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
