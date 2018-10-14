#include <stdio.h>
#include "cblas.h"
#include <lapacke.h>
#include <string.h>

int main (int argc, const char * argv[]) {

  double A[2][2] = {3,4,5,9};
  double B[2][1] = {5,3};

  lapack_int info,m,n=2,lda,ldb,nrhs;
  lapack_int *ipiv;
  int i,j,k;
  
  printf("Input:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     {
        printf("%lf ",A[i][j]);
     }
     printf("  %lf ", B[i][0]);
     printf("\n");
  }
 
  ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  m = n;
  nrhs = 2;
  lda = m;
  ldb = 2;

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, *A, lda, ipiv);
  
  printf("LU:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     {
        printf("%lf ",A[i][j]);
     }
     printf("\n");
  }
  
  double *temp=malloc(sizeof(double));
  for (i=0;i<n;i++)
  {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      memcpy(temp,B[i],sizeof(double));
      memcpy(B[i],B[ipiv[i]],sizeof(double));
      memcpy(B[ipiv[i]],temp,sizeof(double));
    }
  }
  
  cblas_dtrsm(LAPACK_ROW_MAJOR, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, n, 1, *A, n, *B, n);
  printf("y:\n");
  for(i=0;i<n;i++)
  {
     printf("  %lf ", B[i][0]);
     printf("\n");
  }
  cblas_dtrsm(LAPACK_ROW_MAJOR, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1, *A, n, *B, n);
  
  printf("Result:\n");
  for(i=0;i<n;i++)
  {
     printf("  %lf ", B[i][0]);
     printf("\n");
  }
  return 0;
}
