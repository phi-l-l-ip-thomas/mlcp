//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// CUDA-driven "intertwining" solver
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains the CUDA version of the intertwining code for computing eigenvectors

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "utils.h"
#include "matrix.h"
#include "cpvector.h"
#include "hoplist.h"
#include "modvecvec.h"
#include "als.h"
#include "intertw_mvp.h"
#include "intertw_lcv.h"

using namespace std;

void testfltdb(const prec_typ valin);

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

extern"C" void intwdevicemem_(const int &rk, const int &ndof, const int &nbloc,
                                const int &TxD, const int &nops, const int &msize,
                                const int &ncpu, const bool &update, const int *nbas)

////////////////////////////////////////////////////////////////////////
// Prints device memory required by GPU code

{
   int i,nmax,nD;
   int T=TxD/ndof;
   prec_typ MB=1.0/(2 << 19);
   prec_typ HLMB,vecMB,QMB,matMB,ivMB,mvMB,gsMB,geMB,upMB;

// Compute nmax and nD
   for (i=0,nmax=0,nD=0;i<ndof;i++){
       nD+=nbas[i];
       if (nbas[i] > nmax) nmax=nbas[i];
   }

// Device storage for Hamiltonian
   HLMB=((ndof*(T+1)+nops)*sizeof(int) + (T+msize)*sizeof(prec_typ))*MB;

// Device storage for a CP-vector, block of vectors
   vecMB=((2*ndof+1)*sizeof(int) + rk*(nD+1)*sizeof(prec_typ))*MB;
   QMB=nbloc*vecMB;

// Device storage for generalized eigenvalue problem matrices
   matMB=nbloc*nbloc*sizeof(prec_typ)*MB;

// Device storage for common components of 'intvars' structure
   ivMB=((rk*(nmax+rk+1)+1)*sizeof(prec_typ)+(rk+1)*sizeof(int))*MB;

// Device storage for specific components of 'intvars' structure
   mvMB=ncpu*(2*vecMB+ivMB) + HLMB;
   gsMB=(rk*rk+nbloc)*sizeof(prec_typ)*MB + ivMB + QMB;

// Print the available device memory
   printf(" *** Current device memory usage *** \n \n");
   PrintGPUmem();
   printf("\n");

// Print the requirement for each step
   printf(" *** Device memory requirements  *** \n \n");
#ifdef USE_DOUBLES
   printf(" (double precision run) \n");
#else
   printf(" (single precision run) \n");
#endif
   printf(" Hamiltonian storage     : %10.2f MB \n",HLMB);
   printf(" Vector block storage    : %10.2f MB \n",QMB);
   printf(" Mat-vec product   TOTAL : %10.2f MB \n",mvMB);
   printf(" Gram-Schmidt      TOTAL : %10.2f MB \n",gsMB);

   if (update) {
      geMB=ncpu*(2*vecMB+ivMB) + HLMB + QMB + 2*matMB + nbloc*sizeof(prec_typ)*MB;
      upMB=ncpu*(vecMB+ivMB) + QMB + matMB;

      printf(" <Q^T|H|Q> + solve TOTAL : %10.2f MB \n",geMB);
      printf(" Vector updates    TOTAL : %10.2f MB \n",upMB);
   }

   printf("\n");
   printf(" (N.B. excludes CUSOLVER work arrays) \n \n");
}

////////////////////////////////////////////////////////////////////////

extern"C" void intertwinecpp_(prec_typ *base, prec_typ *coef, int *nbas,
                              const int &rk, const int &ndof, const int &nbloc,
                              int *opid, int *mptr, prec_typ *Hcoef, prec_typ *mat,
                              const int &TxD, const int &nops, const int &msize,
                              prec_typ *eigv, const int &npow, const int &ncpu,
                              const int &nals, const bool &update,
                              const int &ishift, const prec_typ &Eshift)

////////////////////////////////////////////////////////////////////////
// Wrapper for performing the intertwined power method on a block of 
// vectors in CP-format, using the GPU. This subroutine does
// 1) Set of npow matrix-vector-products
// 2) Gram-Schmidt orthogonalization of the vectors in the block
// (optional, if 'update' is set to true):
// 3) Solve generalized eigenvalue problem (GEP)
// 4) Replace block of vectors with eigenvectors of the GEP

{
   CPvec Q_h,Q_d;
   Matrix QHQ_d;
   HopList HL_h,HL_d;
   int i;
   int T=TxD/ndof;

// Build nbas with the appropriate length for the Hamiltonian
   int *nbasq = new int [ndof];
   for (i=0;i<ndof;i++){
       nbasq[i]=nbas[i]*(nbas[i]+1)/2;
   }

// Point a CPvec struct to the arrays containing Q
   CPfromComponents(Q_h,base,coef,nbas,ndof,nbloc*rk);

// Point a HopList struct to arrays containing the Hamiltonian
   HLfromComponents(HL_h,T,ndof,nops,msize,nbasq,opid,mptr,Hcoef,mat);

// Send HopList to device for matrix-vector products
   NewHopListonDevice(HL_d,T,ndof,nops,msize);
   SendHopList2GPU(HL_h,HL_d);

// --- MATRIX-VECTOR PRODUCTS ---
   MatrixVectorProducts(HL_d,Q_h,ncpu,npow,nbloc,rk,ishift,Eshift);

// Remove Hamiltonian from device to free up memory for Gram-Schmidt
   FlushHopListonDevice(HL_d);

// Put the block of vectors on device
   SendCPvec2GPU(Q_h,Q_d);

// --- GRAM-SCHMIDT ORTHOGONALIZATION ---
   GramSchmidt(Q_h,Q_d,nals,nbloc,rk);

   if (update){

//    Put Hamiltonian back on device and allocate <Q^T|H|Q> matrix
      NewHopListonDevice(HL_d,T,ndof,nops,msize);
      SendHopList2GPU(HL_h,HL_d);
      NewMatrixonDevice(QHQ_d,nbloc,nbloc);

//    --- GENERALIZED EIGENVALUE PROBLEM ---
      GenEigvalProblem(HL_d,Q_h,Q_d,QHQ_d,eigv,ncpu,nals,nbloc,rk);

      FlushHopListonDevice(HL_d);

//    --- VECTOR UPDATES ---
      VectorUpdates(Q_h,Q_d,QHQ_d,ncpu,nals,nbloc,rk);
      FlushMatrixonDevice(QHQ_d);

   } else {

//    If no updates, just overwrite Q_h on host with orthogonalized Q_d on device
      GetCPvecRtermsfromGPU(Q_h,Q_d,0);
   }

// Clean-up
   FlushCPveconDevice(Q_d);
   delete[] nbasq;
   cudaDeviceReset();
}

////////////////////////////////////////////////////////////////////////

extern"C" void alscpp_(prec_typ *coef, prec_typ *base, int *nbas, 
		       const int &ndof, const int &rF, const int &rG, 
		       const int &nals)

////////////////////////////////////////////////////////////////////////
// Wrapper for performing ALS using the GPU.

{
   CPvec Q_h,Q_d;
   int nbloc=(rF+rG)/rF;

// Point a CPvec struct to the arrays containing F,G; copy to device
   CPfromComponents(Q_h,base,coef,nbas,ndof,rF+rG);
   SendCPvec2GPU(Q_h,Q_d);

// ALS
   ITERATE_ALS(Q_h,Q_d,rF,nbloc,nals);

// Clean-up
   FlushCPveconDevice(Q_d);
   cudaDeviceReset();

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
