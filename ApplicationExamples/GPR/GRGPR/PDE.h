#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void pdevarname_(char* MyNameOUT, int* ind);
void pdeauxname_(char* MyNameOUT, int* ind);
void registerinitialdata_(const char* const id_name, int* id_name_len);
void getnumericalsolution_(double* V,double* Q);
void getexactsolution_(double* V,double* pos,double* timeStamp);

void pderefinecriteria_(int* refine_flag,const double* max_luh,const double* min_luh);

void inittecplot_(const int* N_in,const int* M_in);
void elementcalltecplotplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
void finishtecplotplotter_(const int* Myrank);
void initializetecplotplotter_(const double* time);
}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
