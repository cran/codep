/*************************************************************************
 
 (c) 2008-2023 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Geographic Distances**
 
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

#include"geodists.h"

void dist_geo_hvs(double* from, double* to, int* n, int* tri,
                  double* d,
                  double* r) {
  double* lon[2] = {&from[n[0]], &to[n[1]]};
  double LAT[2], half_delta_lat, half_delta_lon, a,
  sin2_half_delta_lon, radiusx2;
  int end, i, j, k;
  radiusx2 = 2.0*(*r);
  end = (*tri) ? n[0]*(n[0]-1)/2 : n[0]*n[1];
  for(i = 0, j = (*tri) ? 1 : 0, k = 0; k < end; k++) {
    LAT[0] = PI*from[i]/180.0;
    LAT[1] = PI*to[j]/180.0;
    half_delta_lat = 0.5*(LAT[1] - LAT[0]);
    half_delta_lon = 0.5*PI*(lon[1][j] - lon[0][i])/180.0;
    sin2_half_delta_lon = sin(half_delta_lon);
    sin2_half_delta_lon *= sin2_half_delta_lon;
    a = sin(half_delta_lat);
    a *= a;
    a += cos(LAT[0])*cos(LAT[1])*sin2_half_delta_lon;
    a = sqrt(a);
    d[k] = (a < 1.0) ? a*radiusx2 : radiusx2;
    j++;
    if(j == n[1]) {
      i++;
      j = (*tri) ? i+1 : 0;
    }
    if(!(k%INTERRUPT_CHECK))
      R_CheckUserInterrupt();
  }
  return;
}

void dist_geo_vif(double* from, double* to, int* n, int* tri,
                  double* d, int* niter,
                  double* a, double* f, int* maxiter, double* tol) {
  double *lon[2] = {&from[n[0]], &to[n[1]]};
  double b = (1.0 - *f) * (*a);
  double a2minb2overb2 = b*b;
  a2minb2overb2 = ((*a)*(*a) - a2minb2overb2)/a2minb2overb2;
  double U[2], L, lambda, lambda_prev, tmp, C, u2, A, B, k1, deltaSigma;
  double cosU[2], sinU[2], cos_lambda, sin_lambda;
  double cos_sigma, cos_2sigma_m, sin_sigma, sigma;
  double cos2_alpha, sin_alpha;
  double sinU0xcosU1, cosU0xsinU1, cosU1xsin_lambda, sinU0xsinU1,
  cosU0xcosU1, cos_2sigma_m2;
  int end, i, j, k;
  end = (*tri) ? n[0]*(n[0]-1)/2 : n[0]*n[1];
  for(i = 0, j = (*tri) ? 1 : 0, k = 0; k < end; k++) {
    U[0] = atan((1.0 - (*f))*tan(PI*from[i]/180.0));
    U[1] = atan((1.0 - (*f))*tan(PI*to[j]/180.0));
    L = PI*(lon[1][j] - lon[0][i])/180.0;
    lambda = L;
    cosU[0] = cos(U[0]);
    cosU[1] = cos(U[1]);
    sinU[0] = sin(U[0]);
    sinU[1] = sin(U[1]);
    sinU0xcosU1 = sinU[0]*cosU[1];
    cosU0xsinU1 = cosU[0]*sinU[1];
    sinU0xsinU1 = sinU[0]*sinU[1];
    cosU0xcosU1 = cosU[0]*cosU[1];
    // lambda iteration loop:
    do {
      niter[k]++;
      // part 1
      cos_lambda = cos(lambda);
      sin_lambda = sin(lambda);
      cosU1xsin_lambda = cosU[1]*sin_lambda;
      sin_sigma = cosU1xsin_lambda*cosU1xsin_lambda;
      tmp = cosU0xsinU1 - sinU0xcosU1*cos_lambda;
      tmp *= tmp;
      sin_sigma += tmp;
      sin_sigma = sqrt(sin_sigma);
      cos_sigma = sinU0xsinU1 + cosU0xcosU1*cos_lambda;
      sigma = atan2(sin_sigma,cos_sigma);
      sin_alpha = cosU[0]*cosU1xsin_lambda/sin_sigma;
      cos2_alpha = 1.0 - sin_alpha*sin_alpha;
      cos_2sigma_m = cos_sigma - 2.0*sinU0xsinU1/cos2_alpha;
      C = (*f)*cos2_alpha*(4.0 + (*f)*(4.0 - 3.0*cos2_alpha))/16.0;
      lambda_prev = lambda;
      cos_2sigma_m2 = cos_2sigma_m*cos_2sigma_m;
      tmp = cos_2sigma_m + C*cos_sigma*(-1.0 + 2.0*cos_2sigma_m2);
      lambda = L + (1.0 - C)*(*f)*sin_alpha*(sigma + C*sin_sigma*tmp);
    } while ((niter[k] <= *maxiter) && (fabs(lambda - lambda_prev) > *tol));
    u2 = cos2_alpha*a2minb2overb2;
    tmp = sqrt(1.0 + u2);
    k1 = (tmp - 1.0)/(tmp + 1.0);
    tmp = k1*k1;
    A = (1.0 + 0.25*tmp)/(1.0 - k1);
    B = k1*(1.0 - 0.375*tmp);
    tmp = cos_2sigma_m2;
    deltaSigma = B*sin_sigma*
      (cos_2sigma_m + 0.25*B*(cos_sigma*(-1.0 + 2.0*cos_2sigma_m2) -
      B*cos_2sigma_m*(-3.0 + 4.0*sin_sigma*sin_sigma)*
      (-3.0 + 4.0*cos_2sigma_m2)/6));
    d[k] = b*A*(sigma-deltaSigma);
    j++;
    if(j == n[1]) {
      i++;
      j = (*tri) ? i+1 : 0;
    }
    if(!(k%INTERRUPT_CHECK))
      R_CheckUserInterrupt();
  }
  return;
}

/* Computes the Euclidean distance between two matrices.
 The only implementation I could find in R takes a
 single matrix and renders a lower triangular
 distance matrix. That function can take two matrices with
 identical numbers of columns and output a rectangular
 distance matrix. (Wrapper needed here)
 */

void dist_Euclid(double* from, double* to, int* n, int* tri, int* m,
                 double* d,
                 int* sq) {
  int end, i, j, k, l, ii, jj;
  double acc;
  end = (*tri) ? n[0]*(n[0]-1)/2 : n[0]*n[1];
  for(i = 0, j = (*tri) ? 1 : 0, k = 0; k < end; k++) {
    for(l = 0, ii = i, jj = j; l < *m; l++,
        ii += n[0], jj += n[1]) {
      acc = from[ii] - to[jj];
      acc *= acc;
      d[k] += acc;
    }
    if(!*sq)
      d[k] = sqrt(d[k]);
    j++;
    if(j == n[1]) {
      i++;
      j = (*tri) ? i+1 : 0;
    }
    if(!(k%INTERRUPT_CHECK))
      R_CheckUserInterrupt();
  }
  return;
}

