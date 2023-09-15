//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// MODULE linalg
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains wrappers for linear algebra subroutines

#include <iostream>
#include <stdio.h>
#include "cusolverDn.h"
#include "global.h"
#include "utils.h"
#include "matrix.h"
#include "linalg.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void solve_geneigval(Matrix &QHQ_d, const Matrix &S_d, prec_typ *eigv_d, const int &n)

////////////////////////////////////////////////////////////////////////
// Solves general eigenvalue problem. Arrays QHQ_d, S_d, and eigv_d must
// be allocated on device before calling

{
// CUSOLVER handles and initialization
   cusolverDnHandle_t gep;
   cusolverEigType_t ITYP=CUSOLVER_EIG_TYPE_1;
   cusolverEigMode_t JOB=CUSOLVER_EIG_MODE_VECTOR;
   cublasFillMode_t UPLO=CUBLAS_FILL_MODE_LOWER;
   cusolverStatus_t stat;
   stat = cusolverDnCreate(&gep);
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnCreate()");
   }
   prec_typ *WORK;
   int *INFO;
   int LWORK;

// Allocate workspace
#ifdef USE_DOUBLES
   stat = cusolverDnDsygvd_bufferSize(gep,ITYP,JOB,UPLO,n,&QHQ_d.mat[0],n,&S_d.mat[0],n,&eigv_d[0],&LWORK);
#else
   stat = cusolverDnSsygvd_bufferSize(gep,ITYP,JOB,UPLO,n,&QHQ_d.mat[0],n,&S_d.mat[0],n,&eigv_d[0],&LWORK);
#endif
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnXsygvd_bufferSize()");
   }
   cudaMalloc(&WORK,LWORK*sizeof(prec_typ));
   cudaMalloc(&INFO,sizeof(int));

// Solve generalized eigenvalue problem
#ifdef USE_DOUBLES
   stat = cusolverDnDsygvd(gep,ITYP,JOB,UPLO,n,&QHQ_d.mat[0],n,&S_d.mat[0],n,&eigv_d[0],&WORK[0],LWORK,INFO);
#else
   stat = cusolverDnSsygvd(gep,ITYP,JOB,UPLO,n,&QHQ_d.mat[0],n,&S_d.mat[0],n,&eigv_d[0],&WORK[0],LWORK,INFO);
#endif
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnXsygvd()");
   }

// Cleanup
   cudaFree(WORK);
   cudaFree(INFO);
   cusolverDnDestroy(gep);
}

////////////////////////////////////////////////////////////////////////

__host__ void solve_linsys(Matrix &BB_d, Matrix &bjk_d, prec_typ *WORK, int *IPV, int *INFO, const int &rk, const int &n)

////////////////////////////////////////////////////////////////////////
// Solves system of linear equations. All arrays must be allocated on 
// device before calling

{

// Set up CUSOLVER
   cusolverDnHandle_t ls;
   cublasOperation_t trans = CUBLAS_OP_N;
   cusolverStatus_t stat;
   stat = cusolverDnCreate(&ls);
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnCreate()");
   }

// Solve linear system BB*cjk = bjk with DGETRF + DGETRS
#ifdef USE_DOUBLES
   stat = cusolverDnDgetrf(ls,rk,rk,&BB_d.mat[0],rk,&WORK[0],&IPV[0],INFO);
#else
   stat = cusolverDnSgetrf(ls,rk,rk,&BB_d.mat[0],rk,&WORK[0],&IPV[0],INFO);
#endif
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnXgetrf()");
   }

#ifdef USE_DOUBLES
   stat = cusolverDnDgetrs(ls,trans,rk,n,&BB_d.mat[0],rk,&IPV[0],&bjk_d.mat[0],rk,INFO);
#else
   stat = cusolverDnSgetrs(ls,trans,rk,n,&BB_d.mat[0],rk,&IPV[0],&bjk_d.mat[0],rk,INFO);
#endif
   if (stat!=CUSOLVER_STATUS_SUCCESS){
      throwERROR("error with cusolverDnXgetrs()");
   }

// Clean up
   cusolverDnDestroy(ls);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
