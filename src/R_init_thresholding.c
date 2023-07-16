/*
* Code to register compiled routines within R through shared objects
* (DLL/shared object files)
*/

#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <R.h>
#include <stdlib.h> /* Definition for NULL */

/* C/C++ routines with .Call function */
extern SEXP R_thresh_analysis(SEXP, SEXP, SEXP, SEXP, SEXP,
                        SEXP, SEXP, SEXP, SEXP, SEXP,
                        SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallMethods[] = {
    {"R_thresh_analysis", (DL_FUNC) &R_thresh_analysis, 15}, 
    {NULL, NULL, 0}
};

/* Register routines for the packages and allocate resources*/
void R_init_thresholding(DllInfo *info){
    R_registerRoutines(info, NULL, CallMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

/* May not need to implement */
void R_unload_thresholding(DllInfo *info){  ;  }