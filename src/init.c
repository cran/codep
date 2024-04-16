/*************************************************************************
 
 (c) 2008-2023 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Registering routines and dynamic symbols**
 
 This file is part of codep
 
 codep is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 codep is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with codep.  If not, see <https://www.gnu.org/licenses/>.
 
 *************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void mcapermute(double*, double*, double*, int*, double*, double*, int*,
                       int*, int*, int*, int*);
extern void dist_geo_hvs(double*, double*, int*, int*, double*, double*);
extern void dist_geo_vif(double*, double*, int*, int*, double*, int*,
                         double*, double*, int*, double*);
extern void dist_Euclid(double*, double*, int*, int*, int*, double*, int*);
extern void LGTr_C(double*, int*, int*, int*);
// extern void scf_spher(double*, double*, double*, int*, int*, double*);
// extern void scf_expon(double*, double*, double*, int*, int*, double*);
// extern void scf_power(double*, double*, double*, int*, int*, double*);
// extern void scf_hyper(double*, double*, double*, int*, int*, double*);
// extern void scf_super(double*, double*, double*, int*, int*, double*);
// extern void mat_center(double*, int*, int*, double*, double*, double*, int*);
// extern void get_center(double*, int*, int*, double*, double*, double*, int*);
// extern void mat_recenter(double*, int*, int*, double*, double*, double*, int*);
// extern void moran(double*, double*, int*, double*, int*, int*, int*, double*);

static const R_CMethodDef CEntries[] = {
  {"mcapermute",    (DL_FUNC) &mcapermute,   11},
  {"dist_geo_hvs",  (DL_FUNC) &dist_geo_hvs,  6},
  {"dist_geo_vif",  (DL_FUNC) &dist_geo_vif, 10},
  {"dist_Euclid",   (DL_FUNC) &dist_Euclid,   7},
  {"LGTr_C",        (DL_FUNC) &LGTr_C,        4},
//  {"scf_spher",     (DL_FUNC) &scf_spher,     6},
//  {"scf_expon",     (DL_FUNC) &scf_expon,     6},
//  {"scf_power",     (DL_FUNC) &scf_power,     6},
//  {"scf_hyper",     (DL_FUNC) &scf_hyper,     6},
//  {"scf_super",     (DL_FUNC) &scf_super,     6},
//  {"mat_center",    (DL_FUNC) &mat_center,    7},
//  {"get_center",    (DL_FUNC) &get_center,    7},
//  {"mat_recenter",  (DL_FUNC) &mat_recenter,  7},
//  {"moran",         (DL_FUNC) &moran,         8},
  {NULL,            NULL,                     0}
};

static const R_CallMethodDef CallEntries[] = {
  {NULL,            NULL,                     0}
};

void R_init_codep(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
