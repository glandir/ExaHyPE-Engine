#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__


extern "C" {

// FORTRAN functions called by C
void initialdata_(const double* x, const double* const t, double* Q);

// only initialdata_ is used, no more.

void initialalfenwave_(const double* x,const double* const t,  double* Q);
void initialblast_(const double* x, const double* const t, double* Q);
void initialorsagtang_(const double* x,const double* const t,  double* Q);
void initialrotor_(const double* x, const double* const t, double* Q);

// Exact solutions in FORTRAN
//void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);

void initialaccretiondisc_(const double* x, const double* const t,  double* Q);
void initialaccretiondisc3d_(const double* x, const double* const t, double* Q);

}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
