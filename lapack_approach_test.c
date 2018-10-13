#include <stdio.h>
#include <lapacke.h>

int main (int argc, const char * argv[]) {
  
  double A[2][2] = {3,4,5,9};
  double B[2][1] = {5,3};
   
  lapack_int n=2,*ipiv;
  

  n = 2;
  ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, *A, n, ipiv);
  
  int i,j;
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     {
        printf("%lf ",A[i][j]);
     }
     printf("  %6i ", ipiv[i]);
     printf("\n");
  }
  return 0;
}
