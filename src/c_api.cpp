#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "backend/c_api/grplasso.cpp"
#include "backend/c_api/grpSVM.cpp"
#include "backend/c_api/grpLR.cpp"
#include "backend/c_api/grpPR.cpp"


static R_NativePrimitiveArgType grplasso_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, STRSXP, INTSXP, INTSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType grpSVM_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType grpLR_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, STRSXP, REALSXP, REALSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType grpPR_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, STRSXP, REALSXP, REALSXP, INTSXP, REALSXP
};
static const R_CMethodDef c_Methods[] = {
  {"grplasso", (DL_FUNC) &grplasso, 15, grplasso_t},
  {"grpSVM", (DL_FUNC) &grpSVM, 15, grpSVM_t},
  {"grpLR", (DL_FUNC) &grpLR, 17, grpLR_t},
  {"grpPR", (DL_FUNC) &grpPR, 17, grpPR_t},
  {NULL, NULL, 0, NULL}
};
void R_init_myLib(DllInfo *info)
{
  R_registerRoutines(info, c_Methods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
