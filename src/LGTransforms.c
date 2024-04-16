/*************************************************************************
 
 (c) 2008-2023 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Legendre and Gallagher 2001 distance transformations**
 
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
 
 C fonctions definitions
 
 *************************************************************************/

#include"LGTransforms.h"

void LGTr_C(double* x, int* nr, int* nc, int* m) {
  double rsum, *csum, sum;
  int i, j, idx;
  switch(*m) {
  case 1:
    for(i = 0; i < *nr; i++) {
      rsum = 0.0;
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        rsum += x[idx] * x[idx];
      rsum = sqrt(rsum);
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        x[idx] /= rsum;
    }
    break;
  case 2:
    csum = (double*)Calloc(*nc,double);
    sum = 0.0;
    for(i = 0, idx = 0; i < *nc; i++) {
      for(j = 0, csum[i] = 0.0; j < *nr; j++, idx++)
        csum[i] += x[idx];
      sum += csum[i];
      csum[i] = sqrt(csum[i]);
    }
    sum = sqrt(sum);
    for(i = 0; i < *nr; i++) {
      rsum = 0.0;
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        rsum += x[idx];
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        x[idx] = sum * x[idx] / (rsum * csum[j]);
    }
    Free(csum);
    break;
  case 3:
    for(i = 0; i < *nr; i++) {
      rsum = 0.0;
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        rsum += x[idx];
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        x[idx] /= rsum;
    }
    break;
  case 4:
    for(i = 0; i < *nr; i++) {
      rsum = 0.0;
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        rsum += x[idx];
      for(j = 0, idx = i; j < *nc; j++, idx += *nr)
        x[idx] = sqrt(x[idx] / rsum);
    }
    break;
  }
  return;
}
