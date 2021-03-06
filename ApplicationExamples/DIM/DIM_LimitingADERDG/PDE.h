#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void registerinitialdata_(const char* const id_name, int* id_name_len);

void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdejacobian_(double* An, const double* const Q, const double* const gradQ, double* nv);
void pdematrixb_(double* Bn, const double* const Q, double* nv);

void pdeintermediateFields_(double* RL, double* LL, double* iRL, const double* const Q, const double* const nv);
void pdeeigenvectors_(double* L, double* R, double* iR, const double* const Q, const double* const nv);
void hllemriemannsolver_(const int* basisSize, const int* normalNonZeroIndex, double* FL, double* FR, const double* const  QL, const double* const  QR, const double* const  QavL, const double* const  QavR);

}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
