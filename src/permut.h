/*******************************************
*
* C++ headers
*
*******************************************/
#ifndef __PERMUT___HEADERS__
#define __PERMUT___HEADERS__

SEXP transposeMat(SEXP x);
SEXP prodMat(SEXP x, SEXP y);
SEXP prodScal(SEXP x, SEXP y);
SEXP sampleC (SEXP x);
SEXP cbindC(SEXP x, SEXP y);
SEXP selectCol (SEXP x, int y);
SEXP minusMat(SEXP x, SEXP y);
SEXP permut(SEXP p, SEXP ryx, SEXP us, SEXP t2, SEXP details, SEXP permute, SEXP step, SEXP ddfr);

#endif

