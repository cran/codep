/*******************************************************************
 C code to perform permutation testing Multi-scale Codependence
 Analysis (MCA). Handles both univeriate and multivariate testing.
 Guillaume Guenard - Universite de Montreal - 2008-2015
 C functions
*******************************************************************/

// Includes
#include<R.h>
#include<Rmath.h>
#include"codep.h"

// C functions definition
void mcapermute(double *tau0, double *ry, double *rx, double *us, int *n, int *details, int *nperm)
{
  double uspy, uspx, ssqry, ssqrx;
  double rnb, buffer;
  int p, i, j;
  GetRNGstate();                     // This call to insure the RNG is initialized.
  for(p = 0; p < *nperm; p++)
  {
    // 1. Shuffling elements of ry and rx independently
    for(i = 0; i < *n; i++)
    {
      do rnb = unif_rand();
      while (rnb == 1.0);
      j = (int)(rnb * *n);
      buffer = ry[i];
      ry[i] = ry[j];
      ry[j] = buffer;
      do rnb = unif_rand();
      while (rnb == 1.0);
      j = (int)(rnb * *n);
      buffer = rx[i];
      rx[i] = rx[j];
      rx[j] = buffer;
    }
    // 2. Calculation of uspy and uspx
    uspy = 0.0, uspx = 0.0;
    for(i = 0; i < *n; i++)
    {
      uspy += us[i] * ry[i];
      uspx += us[i] * rx[i];
    }
    // 3. Calculation of the residual sum of squares
    ssqry = 0.0, ssqrx = 0.0;
    for(i = 0; i < *n; i++)
    {
      buffer = ry[i] - us[i] * uspy;
      ssqry += buffer * buffer;
      buffer = rx[i] - us[i] * uspx;
      ssqrx += buffer * buffer;
    }
    // 4. Calculation of the permuted tau0
    buffer = uspy * uspx * R_pow(ssqry * ssqrx, -0.5);
    // 5. Comparisons with the nominal value
    *tau0 = ((*tau0 < 0) ? -(*tau0) : *tau0);
    if(buffer <= -(*tau0))
      details[0]++;
    else if(buffer >= *tau0)
      details[2]++;
    else
      details[1]++;
  }
  return;
}

void mmcapermute(double *phi0, double *rY, int *m, double *rx, double *us, int *n, int *details, int *nperm)
{
  double *uspY, uspx, ssqhY, ssqrY, ssqhx, ssqrx;
  double rnb, buffer;
  int p, i, j, k, os1, os2;            // Here, i : rows and j : cols of Y, os1 os2 for table offsetting.
  uspY = (double*)Calloc(*m, double);  // For the 
  GetRNGstate();                       // This call to insure the RNG is initialized.
  for(p = 0; p < *nperm; p++)
  {
  // 1. Shuffling elements of rY and rx independently
    for(i = 0; i < *n; i++)
    {
      do rnb = unif_rand();
      while (rnb == 1.0);
      j = (int)(rnb * *n);
      // Shuffling the rows of rY together.
      for(k = 0, os1 = i, os2 = j; k < *m; k++, os1 += *m, os2 += *m)
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
  // 3. Calcul des sommes de carré....
    ssqhY = 0.0, ssqrY = 0.0;
    // Here, i : rows and j : cols of Y
    for(j = 0, os1 = 0; j < *m; j++)
    {
      for(i = 0; i < *n; i++, os1++)
      {
	buffer = us[i] * uspY[j];      // Calculation of the Yhat[i,j]
	ssqhY += buffer * buffer;      // Accumulation of Yhat's sum of squares
	buffer = rY[os1] - buffer;     // Calculation of residual for Y[i,j]
	ssqrY += buffer * buffer;      // Accumulation of Yres's sum of squares
      }
    }
    ssqhx = 0.0, ssqrx = 0.0;
    for(i = 0; i < *n; i++)
    {
      buffer = us[i] * uspx;
      ssqhx += buffer * buffer;
      buffer = rx[i] - buffer;
      ssqrx += buffer * buffer;
    }
    // 4. Calculation of the permuted tau0
    buffer = ssqhY * ssqhx / (ssqrY * ssqrx);
    // 5. Comparisons with the nominal value
    if(buffer >= *phi0)
      details[1]++;
    else
      details[0]++;
  }  // End of the permutation loop.
  // Free block
  Free(uspY);
  return;
}
