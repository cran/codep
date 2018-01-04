/*******************************************************************
 C code to perform permutation testing Multi-scale Codependence
 Analysis (MCA). Handles both univeriate and multivariate testing.
 Guillaume Guenard - Universite de Montreal - 2008-2018
 C functions
*******************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_codep(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
