//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CUDA-driven "intertwining" solver, Gram-Schmidt/Updates portion
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains the CUDA version of the intertwining code,
// linear-combination-of-vectors part

#include <iostream>
#include <stdio.h>
#include <omp.h>
#include "global.h"
#include "utils.h"
#include "matrix.h"
#include "cpvector.h"
#include "modvecvec.h"
#include "linalg.h"
#include "als.h"
#include "intertw_lcv.h"

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

__global__ void setarray_kernel(const prec_typ E, prec_typ *v, const int n)

////////////////////////////////////////////////////////////////////////
// Sets array 'v' to value 'E'

{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i < n) v[i]=E;
}

////////////////////////////////////////////////////////////////////////

__global__ void gramsumoverrow_kernel(prec_typ *u, const Matrix M1, const Matrix M2, const CPvec v, const CPvec w)

////////////////////////////////////////////////////////////////////////
// Sums elements of PS(!=k)[i,j]*PS(k)[i,j] over j, to compute the GS
// coefs. The dimensions of M1 and M2 must be equal, as this kernel does
// not check.

{
   int nr = M1.nr;
   int row = blockIdx.x * blockDim.x + threadIdx.x;
   int col;

   if (row < nr){
      for (col=0,u[row]=0.0;col<M1.nc;col++){
          u[row]-=w.coef[col]*M1.mat[row+col*nr]*M2.mat[row+col*nr];
      }
      u[row]*=v.coef[row];
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void lcv_build_linear_system(const PVVkernel &psk, const PVVkernel &pst, const PVVkernel &bjk, CPvec &F_d, CPvec &Q_d, Matrix &BB_d, Matrix &BBk_d, Matrix &bjk_d, prec_typ *tmp_d, prec_typ *coef, const int &idof, const int &ni, const int &iv, const bool &gs)

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
   dim3 dB(1);
   dim3 dG(1);
   dim3 dB1(BLOCK_SIZE_1D);
   dim3 dB2(BLOCK_SIZE_2D,BLOCK_SIZE_2D);
   dim3 dG1((rk+dB1.x-1)/dB1.x);                           // rk x 1
   dim3 dimGridbjk((rk+dB2.x-1)/dB2.x,(ni+dB2.y-1)/dB2.y); // rk x n

// Initialize bjk on device
   fillmatrix_kernel<<<dimGridbjk,dB2>>>(bjk_d,beta);

   for (k=0;k<iv;k++){ // Loop over previous vectors
//     Designate the k-th vector in Q as "G_d"
       PointCPvecRtermsonGPU(Q_d,G_d,k*rk,rk);

//     Construct PS(!=idof) matrix
       constpstot_wrapper(pst,G_d,F_d,BB_d,idof,false,false);

       if (gs) {
//        Since Gram-Schmidt needs the full D-dimensional inner product
//        to update the coefs, compute PS(idof) here
          constpsk_wrapper(psk,G_d,F_d,BBk_d,idof,true);

//        Compute the elementwise product of PS=PS(!=idof)*PS(idof) matrices
//        and use the result to update the k-th Gram-Schmidt coefficient
          gramsumoverrow_kernel<<<dG1,dB1>>>(tmp_d,BB_d,BBk_d,G_d,F_d);
          sumarray_kernel<<<dG,dB>>>(&coef[k],tmp_d,rk);
       }

//     Add terms in PS(!=idof) to bjk matrix, weighted by the coef
       constbjk_wrapper(bjk,G_d,BB_d,bjk_d,idof,&coef[k]);
   }

// Build BB(!=idof), which contributes to bjk for GS but not for Vector Updates
   constpstot_wrapper(pst,F_d,F_d,BB_d,idof,false,false);
   if (gs) constbjk_wrapper(bjk,F_d,BB_d,bjk_d,idof,&coef[iv]);
}

////////////////////////////////////////////////////////////////////////

__host__ void GramSchmidt(const CPvec &Q_h, CPvec &Q_d, const int &nals, const int &nbloc, const int &rk)

////////////////////////////////////////////////////////////////////////
// Intertwined Gram-Schmdt orthogonalization. Given CP-format block Q,
// each Q(iv) is iteratively orthogonalized w.r.t all previous Q(iv)

{
   if (nals == 0) return;

   int iv,i,j,n;
   int ndof = Q_d.ndof;
   CPvec F_d;
   intwvars IV;

// Initialize Gram-Schmidt container
   IV.newgs(Q_h,rk,nbloc);

// Initialize CUDA types
   dim3 dimGridb((nbloc+IV.dB1.x-1)/IV.dB1.x); // nbloc x 1

// Set Gram-Schmidt coefs to 1
   setarray_kernel<<<dimGridb,IV.dB1>>>(1.0,IV.coef,nbloc);

// Designate the first vector in Q as "F_d"
   PointCPvecRtermsonGPU(Q_d,F_d,0,rk);

// Normalize the first vector, base and coefs
   constpstot_wrapper(IV.pst,F_d,F_d,IV.BB_d,-1,true,true);// Compute <F_d,F_d>, store in BB_d
   sumoverrow_kernel<<<IV.dG1,IV.dB1>>>(IV.tmp_d,IV.BB_d); // Sum over rows of BB_d
   sumarray_kernel<<<IV.dG,IV.dB>>>(IV.sum_d,IV.tmp_d,rk); // Accumulate the sum over a column
   normalizecoef_kernel<<<IV.dG1,IV.dB1>>>(F_d,IV.sum_d);  // Normalize: divide by sqrt(sum_d)

// Loop over vectors in the block
   for (iv=1;iv<nbloc;iv++){

//     Designate the iv-th vector in Q as "F_d"
       PointCPvecRtermsonGPU(Q_d,F_d,iv*rk,rk);

//     Compute the regularization penalty from the coefs of F
       setpenalty_kernel<<<IV.dG,IV.dB>>>(IV.pen_d,F_d.coef,rk);

//     Loop over intertwining interations
       for (i=0;i<nals;i++){

//         Loop over DOF
           for (j=0;j<ndof;j++){
               n = IV.nbas[j];

//             Calculate BB and b_j_k
               lcv_build_linear_system(IV.psk,IV.pst,IV.bjk,F_d,Q_d,IV.BB_d,IV.BBk_d,IV.bjk_d,IV.tmp_d,IV.coef,j,n,iv,true);

//             Add penalty to BB to avoid ill-conditioning
               addpenalty_kernel<<<IV.dG1,IV.dB1>>>(IV.BB_d,IV.pen_d);

//             Solve linear system B*c_j_k = b_j_k
               solve_linsys(IV.BB_d,IV.bjk_d,IV.WORK,IV.IPV,IV.INFO,rk,n);

//             Construct improved Q(iv) and normalize base
               dim3 dimGridbjk((rk+IV.dB2.x-1)/IV.dB2.x,(n+IV.dB2.y-1)/IV.dB2.y); // rF x n
               replacevwithsolns_kernel<<<dimGridbjk,IV.dB2>>>(F_d,IV.bjk_d,j);
               normonebase_kernel<<<IV.dG1,IV.dB1>>>(F_d,j);

//             Normalize coefficients
               constpstot_wrapper(IV.pst,F_d,F_d,IV.BB_d,-1,true,true);// Reconstruct BB_d
               sumoverrow_kernel<<<IV.dG1,IV.dB1>>>(IV.tmp_d,IV.BB_d); // Sum over rows of BB_d
               sumarray_kernel<<<IV.dG,IV.dB>>>(IV.sum_d,IV.tmp_d,rk); // Accumulate sum over column
               normalizecoef_kernel<<<IV.dG1,IV.dB1>>>(F_d,IV.sum_d);  // Normalize: divide by sqrt(sum_d)
           } // Loop over DOF j
       } // Loop over iterations i
   } // Loop over vectors iv

// Dispose Gram-Schmidt container
   IV.flush();

}

////////////////////////////////////////////////////////////////////////

__host__ void UPDATES_ALS(CPvec &Q_d, intwvars &IV, const int &nals)

////////////////////////////////////////////////////////////////////////
// Intertwined vector updates. Given (old) CP-format block Q, each new
// vector is computed as a rank-reduced linear combination of eigenvectors

{
   if (nals == 0) return;

   int i,j,n;
   int rk = IV.rk;
   int nbloc = IV.nbloc;
   int ndof = Q_d.ndof;
   int nD = IV.nD;

// Initialize CUDA types
   dim3 dimGridb((nbloc+IV.dB1.x-1)/IV.dB1.x); // nbloc x 1
   dim3 dimGridCopy((rk+IV.dB2.x-1)/IV.dB2.x,(nD+IV.dB2.y-1)/IV.dB2.y); // rk x nD

// Copy Fo -> F on device
   copycpcoef_kernel<<<IV.dG1,IV.dB1>>>(IV.Fo_d,IV.F_d);
   copycpbase_kernel<<<dimGridCopy,IV.dB2>>>(IV.Fo_d,IV.F_d);

// Compute the regularization penalty from the coefs of F
   setpenalty_kernel<<<IV.dG,IV.dB>>>(IV.pen_d,IV.F_d.coef,rk);

// Loop over intertwining interations
   for (i=0;i<nals;i++){

//     Loop over DOF
       for (j=0;j<ndof;j++){
           n = IV.nbas[j];

//         Calculate BB and b_j_k
           lcv_build_linear_system(IV.psk,IV.pst,IV.bjk,IV.F_d,Q_d,IV.BB_d,IV.BB_d,IV.bjk_d,IV.tmp_d,IV.coef,j,n,nbloc,false);

//         Add penalty to BB to avoid ill-conditioning
           addpenalty_kernel<<<IV.dG1,IV.dB1>>>(IV.BB_d,IV.pen_d);

//         Solve linear system B*c_j_k = b_j_k
           solve_linsys(IV.BB_d,IV.bjk_d,IV.WORK,IV.IPV,IV.INFO,rk,n);

//         Construct improved Q(ivec) and normalize base
           dim3 dimGridbjk((rk+IV.dB2.x-1)/IV.dB2.x,(n+IV.dB2.y-1)/IV.dB2.y); // rF x n
           replacevwithsolns_kernel<<<dimGridbjk,IV.dB2>>>(IV.F_d,IV.bjk_d,j);
           normonebase_kernel<<<IV.dG1,IV.dB1>>>(IV.F_d,j);
       } // Loop over DOF j
   } // Loop over iterations i

// Normalize coefficients at the end
   constpstot_wrapper(IV.pst,IV.F_d,IV.F_d,IV.BB_d,-1,true,true); // Reconstruct BB_d
   sumoverrow_kernel<<<IV.dG1,IV.dB1>>>(IV.tmp_d,IV.BB_d);   // Sum over rows of BB_d
   sumarray_kernel<<<IV.dG,IV.dB>>>(IV.sum_d,IV.tmp_d,rk);   // Accumulate sum over column
   normalizecoef_kernel<<<IV.dG1,IV.dB1>>>(IV.F_d,IV.sum_d); // Normalize: divide by sqrt(sum_d)

//   gpuErrchk( cudaPeekAtLastError() );   // TEST
//   gpuErrchk( cudaDeviceSynchronize() ); // TEST
}

////////////////////////////////////////////////////////////////////////

__host__ void VectorUpdates(CPvec &Q_h, CPvec &Q_d, const Matrix &QHQ_d, const int &ncpu, const int &nals, const int &nbloc, const int &rk)

////////////////////////////////////////////////////////////////////////
// Intertwined vector updates. Given (old) CP-format block Q, each new
// vector is computed as a rank-reduced linear combination of eigenvectors

{
   int i,j;

// Initialize Updates containers
   intwvars IV[ncpu];
   for (j=0;j<ncpu;j++){
       IV[j].newupdates(Q_h,rk,nbloc);
   }

// Compute new vectors via intertwined updates
#pragma omp parallel shared(IV) private(i,j)
{
#pragma omp for
   for (i=0;i<nbloc;i++){
       j=omp_get_thread_num();
//     Point to a vector on device and corresponding column of the QHQ matrix,
//     compute the updated vector, and then copy new vector to host, overwriting the old vector
       IV[j].point(Q_d,i);
       IV[j].coef=&QHQ_d.mat[i*IV[j].nbloc];
       UPDATES_ALS(Q_d,IV[j],nals);
       IV[j].get(Q_h,i);
   }
}

// Dispose Updates containers
   for (j=0;j<ncpu;j++){
       IV[j].flush();
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
