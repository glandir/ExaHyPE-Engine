#ifndef __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__
#define __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__


extern "C" {

// FORTRAN functions called by C
void initialdata_(const double* x, const double* const t, double* Q);

// only initialdata_ is used, no more.

void InitialPlaneWave_(const double* x, const double* const t, double* Q);
void GaussianBubble_(const double* x, const double* const t, double* Q);

// Smoothing functions for alpha
// void SmoothInterface(const double* alpha, const double* const r, const double* const ICsig);
// Exact solutions in FORTRAN
//void alfenwave_(const double* x, double* Q, const double* /* scalar */ t);


}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */
