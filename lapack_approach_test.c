#include <stdio.h>
#include <lapacke.h>

int main (int argc, const char * argv[]) {
  
   double A[2][2] = {3,4,5,9};
   double B[2][1] = {5,3};
   
   lapack_int info,m,n,lda,ldb,nrhs;
   int i,j;

   m = 2;
   n = 2;
   nrhs = 2;
   lda = 2;
   ldb = 2;

   dgetrf(LAPACK_ROW_MAJOR, m, n,	A, lda);

   for(i=0;i<n;i++)
   {
      for(j=0;j<n;j++)
      {
         printf("%lf ",A[i][j]);
      }
      printf("\n");
   }
   return 0;
}
