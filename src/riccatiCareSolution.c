#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdbool.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <stdio.h>

// Select function for dgees
int sortbynegreal(const double* real, const double* img) {
  if (*real < 0  )
    return(true);
  else
    return(false);
}

// Solves Riccati CARE, given matrix x = {{A, G}, {G, Q}}
SEXP riccatiCareSolution(SEXP x) {
  int *dim, n, m, info, izero = 0, lwork = -1;
  double *work, tmp;
  
  dim = INTEGER(getAttrib(x, R_DimSymbol));

  n = dim[0];
  int bwork[n];
  m = n/2;
    
  dgees_("V", "S", sortbynegreal, dim, (double *) NULL, dim, &izero,
        (double *) NULL, (double *) NULL, (double *) NULL, dim,
        &tmp, &lwork, bwork, &info);
 
  lwork = (int) tmp;
  work = (double*)malloc( lwork*sizeof(double) );
  
  double *wr, *wi, *matz;
  wr = (double*)malloc( n*sizeof(double) );
  wi = (double*)malloc( n*sizeof(double) );
  matz = (double*)malloc( n*n*sizeof(double) );
  
  dgees_("V", "S", sortbynegreal, dim,  REAL(x), dim, &izero, wr, wi,
        matz, dim, work, &lwork, bwork, &info);

  SEXP ans = PROTECT(allocMatrix(REALSXP, m, m));
  
  double *px = REAL(ans), *py;
  py = (double*)malloc( m*m*sizeof(double) );
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      px[i*m+j] = matz[ i + j*n + m ];
      py[i*m+j] = matz[ i + j*n ];
    }  
  }
  
  int *ipiv;
  ipiv = (int*)malloc( m*sizeof(int) );
  int num = m, nrhs = m, lda = m, ldb = m, info2;
  
  dgesv_( &num, &nrhs, py, &lda, ipiv, px, &ldb, &info2);

  free(matz); free(wr); free(wi); free(ipiv);  free(work);
  
  UNPROTECT(1);
  return ans;
}

// Solves Riccati CARE for the inverse, given x = {{A, G}, {G, Q}}
SEXP riccatiCareSolutionInverse(SEXP x) {
  int *dim, n, m, info, izero = 0, lwork = -1;
  double *work, tmp;
  
  dim = INTEGER(getAttrib(x, R_DimSymbol));
  
  n = dim[0];
  int bwork[n];
  m = n/2;
  
  dgees_("V", "S", sortbynegreal, dim, (double *) NULL, dim, &izero,
         (double *) NULL, (double *) NULL, (double *) NULL, dim,
         &tmp, &lwork, bwork, &info);
  
  lwork = (int) tmp;
  work = (double*)malloc( lwork*sizeof(double) );
  
  double *wr, *wi, *matz;
  wr = (double*)malloc( n*sizeof(double) );
  wi = (double*)malloc( n*sizeof(double) );
  matz = (double*)malloc( n*n*sizeof(double) );
  
  dgees_("V", "S", sortbynegreal, dim, REAL(x), dim, &izero, wr, wi,
         matz, dim, work, &lwork, bwork, &info);
  
  SEXP ans = PROTECT(allocMatrix(REALSXP, m, m));
  
  double *py = REAL(ans);
  double *px = (double*)malloc( m*m*sizeof(double) );
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      px[i*m+j] = matz[ i + j*n + m ];
      py[i*m+j] = matz[ i + j*n ];
    }
  }
  
  int *ipiv;
  ipiv = (int*)malloc( m*sizeof(int) );
  int num = m, nrhs = m, lda = m, ldb = m, info2;
  
  dgesv_( &num, &nrhs, px, &lda, ipiv, py, &ldb, &info2);
  
  free(matz); free(wr); free(wi); free(ipiv); free(work);
  
  UNPROTECT(1);
  return ans;
}

// Same as riccatiCareSolutionInverse, but works with Armadillo memptrs.
void care_inv_C(double x[], double result[], int m) {
  int info, izero = 0, lwork = -1;
  double *work, tmp;
  
  int n = m * 2;
  int dim[2] = {n, n};
  int bwork[n];
  
  dgees_("V", "S", sortbynegreal, dim, (double *) NULL, dim, &izero,
         (double *) NULL, (double *) NULL, (double *) NULL, dim,
         &tmp, &lwork, bwork, &info);
  
  lwork = (int) tmp;
  work = (double*)malloc( lwork*sizeof(double) );
  
  double *wr, *wi, *matz;
  wr = (double*)malloc( n*sizeof(double) );
  wi = (double*)malloc( n*sizeof(double) );
  matz = (double*)malloc( n*n*sizeof(double) );
  
  dgees_("V", "S", sortbynegreal, dim, x, dim, &izero, wr, wi,
         matz, dim, work, &lwork, bwork, &info);
  
  double *py = result;
  double *px = (double*)malloc( m*m*sizeof(double) );
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      px[i*m+j] = matz[ i + j*n + m ];
      py[i*m+j] = matz[ i + j*n ];
    }
  }
  
  int *ipiv;
  ipiv = (int*)malloc( m*sizeof(int) );
  int num = m, nrhs = m, lda = m, ldb = m, info2;
  
  dgesv_( &num, &nrhs, px, &lda, ipiv, py, &ldb, &info2);
  
  free(matz); free(wr); free(wi); free(ipiv); free(work);
}

