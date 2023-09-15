// Matrix structure
struct Matrix
{
   int nr,nc;
   prec_typ *mat;
};

// Global code
__global__ void fillmatrix_kernel(Matrix M, const prec_typ E);
__global__ void scalematrix_kernel(Matrix M, const prec_typ E);
__global__ void copymatrix_kernel(Matrix M1, Matrix M2);

// Host code
__host__ void NewMatrixonHost(Matrix &M, const int &nr, const int &nc);
__host__ void FlushMatrixonHost(Matrix &M);
__host__ void NewMatrixonDevice(Matrix &M, const int &nr, const int &nc);
__host__ void FlushMatrixonDevice(Matrix &M);
__host__ void GetMatrixFromDevice(Matrix &M_h, Matrix &M_d);
__host__ void PrintMatrix(const Matrix &M);

// Device code
__device__ prec_typ GetMatrixelement(const Matrix M, const int row, const int col);
__device__ void SetMatrixelement(Matrix M, const int row, const int col, const prec_typ E);
__device__ void UpdateMatrixelement(Matrix M, const int row, const int col, const prec_typ E);
__device__ void DowndateMatrixelement(Matrix M, const int row, const int col, const prec_typ E);
__device__ void AddtoMatrixelement(Matrix M, const int row, const int col, const prec_typ E);

