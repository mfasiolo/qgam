/* 
 Symbol registration initialization: original provided by Brian Ripley.
 Anything called from R should be registered here. (See also NAMESPACE:1)
*/ 
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP qgam_pmmult2(SEXP b, SEXP c, SEXP bt, SEXP ct, SEXP nthreads);

void qgam_pls_fit1(double *y, double *X, double *w,double *wy, double *E, 
                   double *Es, int *n, int *q, int *rE, double *eta,
                   double *penalty, double *rank_tol, int *nt, int *use_wy);

void qgam_gdi2(double *X,double *E,double *Es,double *rS,double *U1,
          double *sp,double *theta,double *z,double *w,double *wz,double *wf,
          double *Dth,double *Det,double *Det2,double *Dth2,double *Det_th,
          double *Det2_th,double *Det3,double *Det_th2,
          double *Det4, double *Det3_th, double *Det2_th2,
          double *beta,double *b1,double *w1,
          double *D1,double *D2,double *P0,double *P1,double *P2,
          double *ldet, double *ldet1,double *ldet2,double *rV,
          double *rank_tol,int *rank_est,
          int *n,int *q, int *M,int *n_theta, int *Mp,int *Enrow,int *rSncol,int *deriv,
          int *fixed_penalty,int *nt,int *type,double *dVkk);

static const R_CallMethodDef CallMethods[] = {
  {"qgam_pmmult2", (DL_FUNC) &qgam_pmmult2, 5},
  {NULL, NULL, 0}
};

R_CMethodDef CEntries[] = { 
  {"qgam_pls_fit1", (DL_FUNC) &qgam_pls_fit1, 14},
  {"qgam_gdi2",(DL_FUNC) &qgam_gdi2,48},
  {NULL, NULL, 0}
};

void R_init_qgam(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
