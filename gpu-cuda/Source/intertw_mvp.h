// Global code
__global__ void mvp_kernel(const HopList HL, const int iterm, const CPvec v, CPvec w);
__global__ void mvp_shared_kernel(const HopList HL, const int iterm, const CPvec v, CPvec w);
__global__ void mvp_shared_kernelB(const HopList HL, const int iterm, const int tilelen, const CPvec v, CPvec w);
__global__ void mvp_registerized_kernel(const HopList HL, const int iterm, const int tilelen, const CPvec v, CPvec w);

// Host code
__host__ void mvp_wrapper(const PVVkernel &PV, const HopList &HL, const int &iterm, const CPvec &v, CPvec &w);
__host__ void mvp_build_linear_system(const PVVkernel mvp, const PVVkernel pst, const PVVkernel bjk, const HopList &HL_d, const CPvec &Fo_d, CPvec &F_d, CPvec &G_d, Matrix &BB_d, Matrix &bjk_d, const int &idof, const int &ni, const int &ishift, const prec_typ &Eshift, const bool recalcB);
__host__ void POW_ALS(const HopList &HL_d, intwvars &IV, const int &npow, const int &ishift, const prec_typ &Eshift);
__host__ void PRODHV_ALS(const HopList &HL_d, intwvars &IV, const int &npow, const int &ishift, const prec_typ &Eshift);
__host__ void MatrixVectorProducts(const HopList &HL_d, CPvec &Q_h, const int &ncpu, const int &npow, const int &nbloc, const int &rk, const int &ishift, const prec_typ &Eshift);
__host__ void GenEigvalProblem(const HopList &HL_d, const CPvec &Q_h, CPvec &Q_d, Matrix &QHQ_d, prec_typ *eigv_h, const int &ncpu, const int &nals, const int &nbloc, const int &rk);

