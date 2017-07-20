#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

typedef SEXP (*sexpPmmultPtr) (SEXP, SEXP, SEXP, SEXP, SEXP);

typedef void (*voidPLS1Ptr) (double*, double*, double*, double*, double*, double*, 
                             int*, int*, int*, double*, double*, double*, int*, int*);

typedef void (*voidGdi2Ptr) (double *,double *,double *,double *,double *, double *,double *,double *,double *,
                             double *,double *, double *,double *,double *,double *,double *, double *,double *,
                             double *, double *, double *, double *, double *,double *,double *, double *,
                             double *,double *,double *,double *, double *, double *, double *,double *,
                             double *,int *, int *, int *, int *,int *, int *,int *,int *,int *, int *, int *, 
                             int *, double *);

SEXP qgam_pmmult2(SEXP b, SEXP c, SEXP bt, SEXP ct, SEXP nthreads) {
  
  static sexpPmmultPtr fun = NULL;
  
  if (fun==NULL) fun = (sexpPmmultPtr) R_GetCCallable("mgcv", "mgcv_pmmult2");
  
  SEXP a = fun(b, c, bt, ct, nthreads);
  
  return(a);
} /* qgam_pmmult2 */



void qgam_pls_fit1(double *y, double *X, double *w,double *wy, double *E, double *Es, 
                   int *n, int *q, int *rE, double *eta,
                   double *penalty, double *rank_tol, int *nt, int *use_wy) {
  
    static voidPLS1Ptr fun = NULL;
  
    if (fun==NULL) fun = (voidPLS1Ptr) R_GetCCallable("mgcv", "pls_fit1");
    
    fun(y, X, w, wy, E, Es, n, q, rE, eta, penalty, rank_tol, nt, use_wy);
    
} /* qgam_pls_fit1 */
  
  
  
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
          int *fixed_penalty,int *nt,int *type,double *dVkk)  {
    
    static voidGdi2Ptr fun = NULL;
    
    if (fun==NULL) fun = (voidGdi2Ptr) R_GetCCallable("mgcv", "gdi2");
    
    fun(X, E, Es, rS,  U1,
          sp, theta, z, w, wz, wf,
          Dth, Det, Det2, Dth2, Det_th,
          Det2_th, Det3, Det_th2,
          Det4,  Det3_th,  Det2_th2,
          beta, b1, w1,
          D1, D2, P0, P1, P2,
          ldet,  ldet1, ldet2, rV,
          rank_tol, rank_est,
          n, q,  M, n_theta,  Mp, Enrow, rSncol, deriv,
          fixed_penalty, nt, type, dVkk);
    
} /* qgam_gdi2 */
  

  
  
  