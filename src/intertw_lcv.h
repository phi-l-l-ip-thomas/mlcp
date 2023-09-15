// Global code
__global__ void setarray_kernel(const prec_typ E, prec_typ *v, const int n);
__global__ void gramsumoverrow_kernel(prec_typ *u, const Matrix M1, const Matrix M2, const CPvec v, const CPvec w);

// Host code
__host__ void lcv_build_linear_system(const PVVkernel &psk, const PVVkernel &pst, const PVVkernel &bjk, CPvec &F_d, CPvec &Q_d, Matrix &BB_d, Matrix &BBk_d, Matrix &bjk_d, prec_typ *tmp_d, prec_typ *coef, const int &idof, const int &ni, const int &iv, const bool &gs);
__host__ void GramSchmidt(const CPvec &Q_h, CPvec &Q_d, const int &nals, const int &nbloc, const int &rk);
__host__ void UPDATES_ALS(CPvec &Q_d, intwvars &IV, const int &nals);
__host__ void VectorUpdates(CPvec &Q_h, CPvec &Q_d, const Matrix &QHQ_d, const int &ncpu, const int &nals, const int &nbloc, const int &rk);
