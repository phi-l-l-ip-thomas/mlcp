// Host code
__host__ void solve_geneigval(Matrix &QHQ_d, const Matrix &S_d, prec_typ *eigv_d, const int &n);
__host__ void solve_linsys(Matrix &BB_d, Matrix &bjk_d, prec_typ *WORK, int *IPV, int *INFO, const int &rk, const int &n);
