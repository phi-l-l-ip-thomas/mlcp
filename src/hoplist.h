// Hamiltonian operator list
struct HopList
{
   int nrk,ndof,nops,msize;
   int *nbas,*opid,*mptr;
   prec_typ *coef,*mat;
};

__host__ void HLfromComponents(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize, int *nbas, int *opid, int *mptr, prec_typ *coef, prec_typ *mat);
__host__ void NewHopListonHost(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize);
__host__ void FlushHopListonHost(HopList &H);
__host__ void NewHopListonDevice(HopList &H, const int &nrk, const int &ndof, const int &nops, const int &msize);
__host__ void FlushHopListonDevice(HopList &H);
__host__ void SendHopList2GPU(HopList &v, HopList &w);
__host__ void GetHopListFromGPU(HopList &v, HopList &w);
__host__ void PrintHopList(const HopList &H);

