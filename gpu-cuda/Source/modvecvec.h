// Inner product kernel structure
struct PVVkernel
{
   dim3 dB,dG;
   int n,m,rkx,rky,rkz,mem,shared,tilelen;
   void initnew(const int &nn, const int &mm, const int &rrkx, const int &rrky, const int &rrkz);
   void init(const int &nn, const int &mm, const int &rrkx, const int &rrky, const int &rrkz);
   void show();
};

// Global code
__global__ void normbase_kernel(CPvec v);
__global__ void normonebase_kernel(CPvec v, const int idof);
__global__ void normalizecoef_kernel(CPvec v, const prec_typ *E);
__global__ void initps_kernel(CPvec v, CPvec w, Matrix PS, const bool multcoef);
__global__ void multcoef_kernel(CPvec v, CPvec w, Matrix PS);
__global__ void constpsk_kernel(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset);
__global__ void constpsk_shared_kernel(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset);
__global__ void constpsk_shared_kernelA(const CPvec v, const CPvec w, Matrix PS, const int k, const bool reset);
__global__ void constpsk_shared_kernelB(const CPvec v, const CPvec w, Matrix PS, const int k, const int tilelen, const bool reset);
__global__ void constpsk_registerized_kernel(const CPvec v, const CPvec w, Matrix PS, const int k, const int tilelen, const bool reset);
__global__ void constpstot_kernel(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef);
__global__ void constpstot_shared_kernel(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef);
__global__ void constpstot_shared_kernelA(CPvec v, CPvec w, Matrix PS, const int idof, const bool multvcoef, const bool multwcoef);
__global__ void constpstot_shared_kernelB(CPvec v, CPvec w, Matrix PS, const int idof, const int tilelen, const bool multvcoef, const bool multwcoef);
__global__ void constpstot_registerized_kernel(CPvec v, CPvec w, Matrix PS, const int idof, const int tilelen, const bool multvcoef, const bool multwcoef);
__global__ void sumoverrow_kernel(prec_typ *v, const Matrix M);
__global__ void sumovercol_kernel(prec_typ *v, const Matrix M);
__global__ void sumarray_kernel(prec_typ *sum, const prec_typ *v, const int n);

// Host code
__host__ void constpsk_wrapper(const PVVkernel &PV, const CPvec &v, const CPvec &w, Matrix &PS, const int &idof, const bool &reset);
__host__ void constpstot_wrapper(const PVVkernel &PV, const CPvec &v, const CPvec &w, Matrix &PS, const int &idof, const bool &multvcoef, const bool &multwcoef);

