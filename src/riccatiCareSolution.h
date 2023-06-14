#ifdef __cplusplus
extern "C" {
#endif
  
  int sortbynegreal(const double* real, const double* img);
  SEXP riccatiCareSolution(SEXP x);
  SEXP riccatiCareSolutionInverse(SEXP x);
  void care_inv_C(double x[], double result[], int m);
  
#ifdef __cplusplus
}
#endif