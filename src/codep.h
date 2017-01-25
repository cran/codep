/*******************************************************************
 C code to perform permutation testing Multi-scale Codependence
 Analysis (MCA). Handles both univeriate and multivariate testing.
 Guillaume Guenard - Universite de Montreal - 2008-2015
 C header
*******************************************************************/

// Defines (empty)

// Type declarations (empty)

// C functions declaration
void mcapermute(double *phi_global0, double *tau_ind0, double *rY, int *m, double *rx, double *us, int *n, int *perm_global, int *perm_ind, int *nperm, int *ind);
