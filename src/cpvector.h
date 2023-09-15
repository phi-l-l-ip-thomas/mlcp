// CP-vector structure
struct CPvec
{
   int nrk,ndof,nD,nmax;
   int *nbas,*nrdim;
   prec_typ *coef;
   prec_typ *base;
};

// Global code
__global__ void copycpcoef_kernel(const CPvec v, CPvec w);
__global__ void copycpbase_kernel(const CPvec v, CPvec w);
__global__ void setcpcoef_kernel(CPvec v, const prec_typ E);
__global__ void scalecpcoef_kernel(CPvec v, const prec_typ E);
__global__ void scalecpbase_kernel(CPvec v, const int idof, const prec_typ E);

// Host code
__host__ void CPfromComponents(CPvec &v, prec_typ *base, prec_typ *coef, int *nbas, const int &ndof, const int &nrk);
__host__ void NewCPveconHost(CPvec &v, const int *nbas, const int &ndof, const int &nrk);
__host__ void FlushCPveconHost(CPvec &v);
__host__ void NewCPveconDevice(CPvec &v, const int *nbas, const int &ndof, const int &nrk);
__host__ void FlushCPveconDevice(CPvec &v);
__host__ void SendCPvec2GPU(const CPvec &v, CPvec &w);
__host__ void SendCPvecRterms2GPU(const CPvec &v, CPvec &w, const int &iR, const int &R);
__host__ void PointCPvecRtermsonGPU(CPvec &v, CPvec &w, const int &iR, const int &R);
__host__ void GetCPvecFromGPU(CPvec &v, const CPvec &w);
__host__ void GetCPvecRtermsfromGPU(CPvec &v, const CPvec &w, const int &iR);
__host__ void PrintCPvec(const CPvec &v);
__host__ int getsharedmemCP(const dim3 &dB, const int &n, const int &m);

// Device code
__device__ prec_typ GetCPelement(const CPvec v, const int irk, const int id, const int ibas);
__device__ void SetCPelement(CPvec v, const int irk, const int idof, const int ibas, const prec_typ E);
__device__ void MultiplyCPelement(CPvec v, const int irk, const int idof, const int ibas, const prec_typ E);

