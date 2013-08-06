#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <stddef.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "permut.h"

#ifdef	__cplusplus
extern "C" {
#endif
	
// Computation of t(X)
SEXP transposeMat(SEXP x)
{
	double *rx = REAL(x), *rz;
	SEXP z, dim;
	R_len_t i, j, xCols, xRows;
	
	PROTECT(z = allocVector(REALSXP,length(x)));
	rz = REAL(z);

	xRows = nrows(x);
	xCols = ncols(x);

	for(i = 0; i < xCols; i++) 
	{		
		for(j = 0; j < xRows; j++)
		{
			rz[i + xCols*j] = rx[xRows*i+j];
		}
	}
	PROTECT(dim = allocVector(INTSXP, 2));
	INTEGER(dim)[0] = xCols;
	INTEGER(dim)[1] = xRows;
	setAttrib(z, R_DimSymbol, dim);
	UNPROTECT(2);
	return(z);
}

// Computation of X %*% Y
SEXP prodMat(SEXP x, SEXP y)
{
	
	double *rx = REAL(x);
	R_len_t i, j, k, l, indice, xRows, xCols, yRows, yCols;
	double *ry = REAL(y), *rz, pdt, tmp;
	
	SEXP z, dimt;
	
	xRows = nrows(x);
	xCols = ncols(x);
	yRows = nrows(y);
	yCols = ncols(y);
	
	if (xCols != yRows)
	{
		error("The number of columns of the first matrix does not match the number of rows of the second matrix");
	}
	
	PROTECT(z = allocVector(REALSXP,xRows*yCols));
	rz = REAL(z);
	indice = 0;

	for (i=0 ; i < yCols && indice < xRows*yCols ; i++)
	{
		for (l=0; l < xRows && indice < xRows*yCols ; l++)
		{
			
			k = 0;
			pdt = 0.0;
			for (j=0 ; j < xCols && indice < xRows*yCols ; j++)
			{				
				tmp = rx[(j*xRows)+l] * ry[k+(yRows*i)];
				pdt = pdt + tmp;
				k++;
			}
			rz[indice] = pdt;
			indice++;
		}
	}
	PROTECT(dimt = allocVector(INTSXP, 2));
	INTEGER(dimt)[0] = xRows; INTEGER(dimt)[1] = yCols;
	setAttrib(z, R_DimSymbol, dimt);
	UNPROTECT(2);
	return(z);
}

// Computation of x * Y
SEXP prodScal(SEXP x, SEXP y)
{
	R_len_t  i;
	double *rx = REAL(x);	
	double *ry = REAL(y), *rz;
	SEXP z, dimps;
	
	
	PROTECT(z = allocVector(REALSXP,length(y)));
	rz = REAL(z);
	
	for (i = 0 ; i < length(y) ; i++)
	{
		rz[i] = rx[0] * ry[i];
	}
	
	PROTECT(dimps = allocVector(INTSXP, 2));
	INTEGER(dimps)[0] = nrows(y); INTEGER(dimps)[1] = ncols(y);
	setAttrib(z, R_DimSymbol, dimps);
	UNPROTECT(2);
	return (z);
}

#define FOR_RAND 1/RAND_MAX
// Computation of sample(x)
SEXP sampleC (SEXP x)
{
	R_len_t i, lx, rdm;
	double *rx = REAL(x), *rz, tmp;
	SEXP z, dims;
	lx = length (x);
	PROTECT(z = allocVector(REALSXP,lx));
	rz = REAL(z);

	for (i = 0 ; i <lx ; i++)
	{
		rz[i] = rx[i];
	}
	for (i = 0 ; i < lx ; i++)
	{
		do {
			rdm = rand();
		} while (rdm == RAND_MAX);
		rdm = (R_len_t) ((double) rdm * (lx) * FOR_RAND);
		tmp = rz[rdm];
		rz[rdm]  = rz[i];
		rz[i] = tmp;
	}
	PROTECT(dims = allocVector(INTSXP, 2));
	INTEGER(dims)[0] = nrows(x); INTEGER(dims)[1] = ncols(x);
	setAttrib(z, R_DimSymbol, dims);
	UNPROTECT(2);
	return(z);
}

// Computation of cbind(x,y)
SEXP cbindC(SEXP x, SEXP y)
{
	R_len_t i;
	SEXP z, dimc;

	PROTECT(z = allocVector(REALSXP,length(x)+length(y)));
	
	for (i = 0 ; i < length(x) + length(y) ; i++)
	{
		if (i < length(x))
		{
			REAL(z)[i] = REAL(x)[i];
		}
		else
		{
			REAL(z)[i] = REAL(y)[i-length(x)];
		}
	}
	PROTECT(dimc = allocVector(INTSXP, 2));
	INTEGER(dimc)[0] = nrows(x); INTEGER(dimc)[1] = ncols(x)+ncols(y);
	setAttrib(z, R_DimSymbol, dimc);
	UNPROTECT(2);
	return(z);
}

// Computation of x[ ,y]
SEXP selectCol (SEXP x, int y)
{
	R_len_t i;
	double *rx = REAL(x), *rz;
	SEXP z, dimsc;
	PROTECT(z = allocVector(REALSXP,nrows(x)));
	rz = REAL(z);

	
	for (i = 0 ; i < nrows(x) ; i++)
	{
		rz[i] = rx[(y-1)*nrows(x) + i];
	}
	
	PROTECT(dimsc = allocVector(INTSXP, 2));
	INTEGER(dimsc)[0] = nrows(x); INTEGER(dimsc)[1] = 1;
	setAttrib(z, R_DimSymbol, dimsc);
	UNPROTECT(2);
	return(z);
}

// Computation of X - Y
SEXP minusMat(SEXP x, SEXP y)
{
	R_len_t i;
	double *rx = REAL(x);
	double *ry = REAL(y), *rz;
	SEXP z, dimmm;
		
	PROTECT(z = allocVector(REALSXP,length(x)));
	rz = REAL(z);

	if (ncols(x) != ncols(y) || nrows(x) != nrows(y))
	{
		error("The first matrix do not match the second matrix");
	}
	rz[0] = 1;
	for (i = 0 ; i < length(x) ; i++)
	{
		rz[i] = rx[i] - ry[i];
	}
	PROTECT(dimmm = allocVector(INTSXP, 2));
	INTEGER(dimmm)[0] = nrows(x); INTEGER(dimmm)[1] = ncols(x);
	setAttrib(z, R_DimSymbol, dimmm);
	UNPROTECT(2);
	return(z);
}

// computation by permutation
SEXP permut(SEXP p, SEXP ryx, SEXP us, SEXP t2, SEXP details, SEXP permute, SEXP step, SEXP ddfr)
{
	R_len_t i;
	SEXP z, p12, p34, uspyxi, t2i, tmp1, tmp2, tmp3, tmp4, tmp5;
	PROTECT_INDEX iz, ip12, ip34, iuspyxi, it2i, itmp1, itmp2, itmp3, itmp4, itmp5;
	SEXP p3, p4, ryx1, ryx2, tus, pmus, uspyxi1, uspyxi2;
	PROTECT_INDEX ip3, ip4, iryx1, iryx2, itus, ipmus , iuspyxi1, iuspyxi2;

	PROTECT_WITH_INDEX(z = allocMatrix(REALSXP,nrows(details),ncols(details)), &iz);
	PROTECT_WITH_INDEX(p12 = allocMatrix(REALSXP,nrows(p),2), &ip12);
	PROTECT_WITH_INDEX(p34 = allocMatrix(REALSXP,nrows(p),2), &ip34);
	PROTECT_WITH_INDEX(uspyxi = allocVector(REALSXP,ncols(us)*2), &iuspyxi);
	PROTECT_WITH_INDEX(t2i = allocVector(REALSXP,1), &it2i);
	PROTECT_WITH_INDEX(tmp1 = allocVector(REALSXP,1), &itmp1);
	PROTECT_WITH_INDEX(tmp2 = allocVector(REALSXP,1), &itmp2);
	PROTECT_WITH_INDEX(tmp3 = allocVector(REALSXP,nrows(details)), &itmp3);
	PROTECT_WITH_INDEX(tmp4 = allocVector(REALSXP,nrows(details)), &itmp4);
	PROTECT_WITH_INDEX(tmp5 = allocVector(REALSXP,nrows(details)), &itmp5);
	PROTECT_WITH_INDEX(p3 = allocVector(REALSXP,1), &ip3);
	PROTECT_WITH_INDEX(p4 = allocVector(REALSXP,1), &ip4);
	PROTECT_WITH_INDEX(ryx1 = allocVector(REALSXP,nrows(ryx)), &iryx1);
	PROTECT_WITH_INDEX(ryx2 = allocVector(REALSXP,nrows(ryx)), &iryx2);
	PROTECT_WITH_INDEX(tus = allocVector(REALSXP,length(us)), &itus);
	PROTECT_WITH_INDEX(pmus = allocVector(REALSXP,length(p)/2), &ipmus);
	PROTECT_WITH_INDEX(uspyxi1 = allocVector(REALSXP,ncols(us)),&iuspyxi1);
	PROTECT_WITH_INDEX(uspyxi2 = allocVector(REALSXP,ncols(us)),&iuspyxi2);
	for (i = 0 ; i < length(details); i++)
	{
		REAL(z)[i] = REAL(details)[i];
	}
	
	for (i = 0 ; i < REAL(permute)[0] ; i++)
	{
		// Computation of p[,1:2]
		REPROTECT(ryx1 = selectCol(ryx,1), iryx1);
		REPROTECT(ryx2 = selectCol(ryx,1), iryx2);
		REPROTECT(ryx1 = sampleC(ryx1), iryx1);
		REPROTECT(ryx2 = sampleC(ryx2), iryx2);
		REPROTECT(p12 = cbindC(ryx1,ryx2), ip12);
		// Computation of uspyxi
		REPROTECT(tus = transposeMat(us), itus);
		REPROTECT(uspyxi = prodMat(tus,p12), iuspyxi);
		// Computation of p[,3:4]
		REPROTECT(pmus = prodMat(us,uspyxi), ipmus);
		REPROTECT(p34 = minusMat(p12,pmus), ip34);
		
		// Computation of t2i
		REPROTECT(uspyxi1 = selectCol(uspyxi,1), iuspyxi1);
		REPROTECT(uspyxi1 = transposeMat(uspyxi1),iuspyxi1);
		REPROTECT(tmp1 = selectCol(uspyxi,REAL(step)[0]), itmp1);
		REPROTECT(uspyxi2 = selectCol(uspyxi,2), iuspyxi2);
		REPROTECT(uspyxi2 = transposeMat(uspyxi2), iuspyxi2);
		REPROTECT(tmp2 = selectCol(uspyxi2,REAL(step)[0]), itmp2);
		REAL(tmp1)[0] = REAL(tmp1)[0]* REAL(tmp2)[0];
		
		REPROTECT(p3 = selectCol(p34,1), ip3);
		REPROTECT(p3 = prodMat(transposeMat(p3),p3), ip3);
		REPROTECT(p4 = selectCol(p34,2), ip4);
		REPROTECT(p4 = prodMat(transposeMat(p4),p4), ip4);
		
		REAL(tmp2)[0] = REAL(p3)[0] * REAL(p4)[0];
		REAL(tmp1)[0] = REAL(tmp1)[0] / sqrt(REAL(tmp2)[0]);
		REAL(t2i)[0] = REAL(ddfr)[0]*REAL(tmp1)[0] ;
		
		// Computation of matrix details
		REPROTECT(tmp3 = selectCol(z,1), itmp3);
		REPROTECT(tmp4 = selectCol(z,2), itmp4);
		REPROTECT(tmp5 = selectCol(z,3), itmp5);
		
		if (REAL(t2i)[0] <= -fabs(REAL(t2)[0]))
		{
			REAL(tmp3)[(int)REAL(step)[0]-1] = REAL(tmp3)[(int)REAL(step)[0]-1] + 1;
		}
		else
		{
			if (REAL(t2i)[0] >= fabs(REAL(t2)[0]))
			{
				REAL(tmp5)[(int)REAL(step)[0]-1] = REAL(tmp5)[(int)REAL(step)[0]-1] + 1;
			}
			else
			{
				REAL(tmp4)[(int)REAL(step)[0]-1] = REAL(tmp4)[(int)REAL(step)[0]-1] + 1;
			}
		}
		REPROTECT(z = cbindC(tmp3,tmp4) , iz);
		REPROTECT(z = cbindC(z,tmp5) , iz);
	}
	

	
	UNPROTECT(18);
	return (z);	
		
}


#ifdef	__cplusplus
}
#endif

