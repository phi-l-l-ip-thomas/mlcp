// Container for intertwining objects
struct intwvars
{
   int nbloc,rk,nD,ndof,nmax,which;
   int *nbas,*INFO,*IPV;
   prec_typ *WORK,*tmp_d,*sum_d,*pen_d,*coef;
   dim3 dB,dB1,dB2,dG,dG1,dG2;
   PVVkernel mvp,psk,pst,bjk;
   Matrix BB_d,BBk_d,bjk_d;
   CPvec Fo_d,F_d,G_d;
   void init(const CPvec &Q_h, const int &rrk);
   void newmvp(const CPvec &Q_h, const int &rrk);
   void newgep(const CPvec &Q_h, const int &rrk, const int &nnbloc);
   void newgs(const CPvec &Q_h, const int &rrk, const int &nnbloc);
   void newupdates(const CPvec &Q_h, const int &rrk, const int &nnbloc);
   void newals(const CPvec &Q_h, const int &rrk, const int &nnbloc);
   void send(const CPvec &Q_h, const int &ivec);
   void get(CPvec &Q_h, const int &ivec);
   void point(CPvec &Q_d, const int &ivec);
   void flush();
};

// Global code
__global__ void setpenalty_kernel(prec_typ *pen, const prec_typ *v, const int n);
__global__ void addpenalty_kernel(Matrix M, const prec_typ *E);
__global__ void addpenalty_kernel(Matrix M, const prec_typ E);
__global__ void replacevwithsolns_kernel(CPvec v, const Matrix bjk, const int idof);
__global__ void constbjk_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k);
__global__ void constbjk_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const prec_typ *fac);
__global__ void constbjk_shared_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k);
__global__ void constbjk_shared_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const prec_typ *fac);
__global__ void constbjk_shared_kernelB(const CPvec w, const Matrix PS, Matrix bjk, const int k, const int tilelen);
__global__ void constbjk_shared_kernelB(const CPvec w, const Matrix PS, Matrix bjk, const int k, const int tilelen, const prec_typ *fac);
__global__ void constbjk_registerized_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const int tilelen);
__global__ void constbjk_registerized_kernel(const CPvec w, const Matrix PS, Matrix bjk, const int k, const int tilelen, const prec_typ *fac);

// Host code
__host__ void als_build_linear_system(const PVVkernel &pst, const PVVkernel &bjk, CPvec &F_d, CPvec &Q_d, Matrix &BB_d, Matrix &bjk_d, const int &idof, const int &ni, const int &iv);
__host__ void ITERATE_ALS(CPvec &Q_h, CPvec &Q_d, const int &rk, const int &nbloc, const int &nals);
__host__ void constbjk_wrapper(const PVVkernel &PV, const CPvec &w, const Matrix &PS, Matrix &bjk, const int &k);
__host__ void constbjk_wrapper(const PVVkernel &PV, const CPvec &w, const Matrix &PS, Matrix &bjk, const int &k, const prec_typ *fac);

