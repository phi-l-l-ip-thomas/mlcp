//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// Matrix operations
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This module contains memory and I/O operations on 2D arrays of class Matrix
// PST, last update 05/30/2017

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "matrix.h"

using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Global functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__global__ void fillmatrix_kernel(Matrix M, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Fills matrix with value 'E', device kernel

{
// ir=row index, ic=col index, of M
   int ir = blockIdx.x * blockDim.x + threadIdx.x;
   int ic = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread writes a value to M
   if (ir < M.nr && ic < M.nc){
      SetMatrixelement(M,ir,ic,E);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void scalematrix_kernel(Matrix M, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Multiplies matrix by value 'E', device kernel

{
// ir=row index, ic=col index, of M
   int ir = blockIdx.x * blockDim.x + threadIdx.x;
   int ic = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread multiplies a value to M by 'E'
   if (ir < M.nr && ic < M.nc){
      UpdateMatrixelement(M,ir,ic,E);
   }
}

////////////////////////////////////////////////////////////////////////

__global__ void copymatrix_kernel(Matrix M1, Matrix M2)

////////////////////////////////////////////////////////////////////////
// Copies elements of M1 -> M2 on device. M1 and M2 should have same
// dimensions as this subroutine does not check

{
// ir=row index, ic=col index, of M1
   int ir = blockIdx.x * blockDim.x + threadIdx.x;
   int ic = blockIdx.y * blockDim.y + threadIdx.y;

// Each thread writes a value to M1
   if (ir < M1.nr && ic < M1.nc){
      int i = ir + ic * M1.nr;
      M2.mat[i] = M1.mat[i];
   }
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Host-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__host__ void NewMatrixonHost(Matrix &M, const int &nr, const int &nc)

////////////////////////////////////////////////////////////////////////
// Initialize Matrix on host

{
   M.nr = nr;
   M.nc = nc;
   M.mat = new prec_typ [nr*nc];
}

////////////////////////////////////////////////////////////////////////

__host__ void FlushMatrixonHost(Matrix &M)

////////////////////////////////////////////////////////////////////////
// Deallocate matrix on host

{
   M.nr=0;
   M.nc=0;
   delete[] M.mat;
}

////////////////////////////////////////////////////////////////////////

__host__ void NewMatrixonDevice(Matrix &M, const int &nr, const int &nc)

////////////////////////////////////////////////////////////////////////
// Initialize Matrix on device

{
   M.nr = nr;
   M.nc = nc;
   cudaMalloc(&M.mat,nr*nc*sizeof(prec_typ));
}

////////////////////////////////////////////////////////////////////////

__host__ void FlushMatrixonDevice(Matrix &M)

////////////////////////////////////////////////////////////////////////
// Deallocate Matrix on device

{
   M.nr = 0;
   M.nc = 0;
   cudaFree(M.mat);
}

////////////////////////////////////////////////////////////////////////

__host__ void GetMatrixFromDevice(Matrix &M_h, Matrix &M_d)

////////////////////////////////////////////////////////////////////////
// Copies Matrix 'M_h' on host <- 'M_d' on device
// 'M_h' should be allocated on host before calling

{
   cudaMemcpy(M_h.mat,M_d.mat,M_h.nr*M_h.nc*sizeof(prec_typ),cudaMemcpyDeviceToHost);
}

////////////////////////////////////////////////////////////////////////

__host__ void PrintMatrix(const Matrix &M)

////////////////////////////////////////////////////////////////////////
// Prints matrix

{
   int i,j;

   for (i=0;i<M.nr;i++){
       for (j=0;j<M.nc;j++){
           printf(" %15.10f ",M.mat[i+j*M.nr]);
       }
       cout << endl;
   }

}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Device-callable functions
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////////////////////////////////////////////////////////////////////////

__device__ prec_typ GetMatrixelement(const Matrix M, const int row, const int col)

////////////////////////////////////////////////////////////////////////
// Retrieves a Matrix element from device memory

{
   return M.mat[row+col*M.nr];
}

////////////////////////////////////////////////////////////////////////

__device__ void SetMatrixelement(Matrix M, const int row, const int col, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Assigns a Matrix element on device memory

{
   M.mat[row+col*M.nr]=E;
}

////////////////////////////////////////////////////////////////////////

__device__ void UpdateMatrixelement(Matrix M, const int row, const int col, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Updates a matrix element my multiplication

{
   M.mat[row+col*M.nr]*=E;
}

////////////////////////////////////////////////////////////////////////

__device__ void DowndateMatrixelement(Matrix M, const int row, const int col, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Downdates a matrix element my division

{
   M.mat[row+col*M.nr]/=E;
}

////////////////////////////////////////////////////////////////////////

__device__ void AddtoMatrixelement(Matrix M, const int row, const int col, const prec_typ E)

////////////////////////////////////////////////////////////////////////
// Adds 'E' to matrix element

{
   M.mat[row+col*M.nr]+=E;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
