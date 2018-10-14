#include <stdio.h>
#include "cblas.h"
#include <lapacke.h>
#include <string.h>

int main (int argc, const char * argv[]) {

  double A[2][2] = {1,0,1,1};
  double B[2][1] = {1,1};

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

  
  cblas_dtrsm(LAPACK_ROW_MAJOR, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, n, 1, *A, n, *B, n);
  printf("Result:\n");
  for(i=0;i<n;i++)
  {
     printf("  %lf ", B[i][0]);
     printf("\n");
  }

  return 0;
}
