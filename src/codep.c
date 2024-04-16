/*************************************************************************
 
 (c) 2008-2023 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Permutation Test for Multiscale Codependence Analysis**
 
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
 
 C functions definitions
 
 *************************************************************************/

// Includes
#include"codep.h"

// C functions definition
void mcapermute(double *phi_global0, double *tau_ind0, double *rY, int *m,
                double *rx, double *us, int *n, int *perm_global, int *perm_ind,
                int *nperm, int *ind)
{
  double *uspY, uspx, *ssqhYi, ssqhY, *ssqrYi, ssqrY, ssqhx, ssqrx;
  double rnb, buffer;
  // Here, i : rows and j : cols of Y; os1, os2, and os3 for table offsetting.
  int p, i, j, k, os1, os2, os3 = *m + *m;
  uspY = (double*)Calloc(*m, double);
  ssqhYi = (double*)Calloc(*m, double);
  ssqrYi = (double*)Calloc(*m, double);
  GetRNGstate();                // This call to insure the RNG is initialized.
  for(p = 0; p < *nperm; p++)
  {
  // 1. Shuffling elements of rY and rx independently
    for(i = 0; i < *n; i++)
    {
      do rnb = unif_rand();
      while (rnb == 1.0);
      j = (int)(rnb * *n);
      // Shuffling the rows of rY together.
      for(k = 0, os1 = i, os2 = j; k < *m; k++, os1 += *n, os2 += *n)
      {
        buffer = rY[os1];
        rY[os1] = rY[os2];
        rY[os2] = buffer;
      }
      do rnb = unif_rand();
      while (rnb == 1.0);
      j = (int)(rnb * *n);
      buffer = rx[i];
      rx[i] = rx[j];
      rx[j] = buffer;
    }
  // 2. Calculation of uspY and uspx.
    // Here, i : rows and j : cols of Y
    for(j = 0, os1 = 0; j < *m; j++)
    {
      uspY[j] = 0.0;
      for(i = 0; i < *n; i++, os1++)
        uspY[j] += us[i] * rY[os1];
    }
    uspx = 0.0;
    for(i = 0; i < *n; i++)
      uspx += us[i] * rx[i];
  // 3. Calculation of sums of squares.
    ssqhY = 0.0, ssqrY = 0.0;
    // Here, i : rows and j : cols of Y
    for(j = 0, os1 = 0; j < *m; j++)
    {
      ssqhYi[j] = 0.0, ssqrYi[j] = 0.0;
      for(i = 0; i < *n; i++, os1++)
      {
	buffer = us[i] * uspY[j];        // Calculation of the Yhat[i,j]
	ssqhYi[j] += buffer * buffer;    // Accumulation of Yhat's sum of squares
	buffer = rY[os1] - buffer;       // Calculation of residual for Y[i,j]
	ssqrYi[j] += buffer * buffer;    // Accumulation of Yres's sum of squares
      }
      ssqhY += ssqhYi[j], ssqrY += ssqrYi[j];
    }
    ssqhx = 0.0, ssqrx = 0.0;
    for(i = 0; i < *n; i++)
    {
      buffer = us[i] * uspx;
      ssqhx += buffer * buffer;
      buffer = rx[i] - buffer;
      ssqrx += buffer * buffer;
    }
  // 4. Calculation of the permuted phi_global0* et phi_resp0*
    buffer = ssqhY * ssqhx / (ssqrY * ssqrx);
  // 5. Comparisons with the nominal value
    // for the global test
    if(buffer >= *phi_global0)
      perm_global[1]++;
    else
      perm_global[0]++;
    // for the individual tests, if requested.
    if(*ind)
      for(j = 0; j < *m; j++)
      {
        buffer = uspY[j] * uspx * R_pow(ssqrYi[j] * ssqrx, -0.5);
	if(buffer <= -tau_ind0[j])
	  perm_ind[j]++;
        else if(buffer >= tau_ind0[j])
	  perm_ind[j + os3]++;
        else
	  perm_ind[j + *m]++;
      }
  }  // End of the permutation loop.
  // Free block
  Free(ssqrYi);
  Free(ssqhYi);
  Free(uspY);
  return;
}
