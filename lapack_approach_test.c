#include <stdio.h>
#include "cblas.h"
#include <lapacke.h>
#include <string.h>

int main (int argc, const char * argv[]) {

  double A[2][2] = {3,4,5,9};
  double B[2][1] = {5,3};

  lapack_int info,m,n,lda,ldb,nrhs;
  lapack_int *ipiv;
  int i,j;

  n = 2;
  ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  m = n;
  nrhs = 2;
  lda = m;
  ldb = 2;

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, *A, lda, ipiv);
  double *temp=malloc(n);

  for (i=0;i<n;i++)
  {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      memcpy(temp,A[i],n);
      memcpy(A[i],A[ipiv[i]],n);
      memcpy(A[ipiv[i]],temp,n);
    }
  }

  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      printf("%lf ",A[i][j]);
    }
    printf("  %lf ", B[i][0]);
    printf("\n");
  }
  
  cblas_dtrsm(LAPACK_ROW_MAJOR, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, *A, n, *B, n);

  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     {
        printf("%lf ",A[i][j]);
     }
     printf("  %lf ", B[i][0]);
     printf("\n");
  }
  return 0;
}
