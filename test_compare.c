#include <stdio.h>
#include <lapacke.h>
#include <time.h>
#include <stdbool.h>
#include <cblas.h>
#include <math.h>
#include <string.h>

int main (int argc, const char * argv[]) {

  int n=3;
  bool test=true;

  double *A=(double *)malloc(sizeof(double)*n*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  int i, j, k;
  A[0]=1;
  A[1]=1;
  A[2]=1;
  A[3]=2;
  A[4]=3;
  A[5]=1;
  A[6]=4;
  A[7]=3;
  A[8]=2;
  B[0]=6;
  B[1]=11;
  B[2]=16;
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind;
  double sum, max;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  clock_t start=clock();
  int *pvt = (int *)malloc(sizeof(int)*n);
  for (i=0;i<n;i++)
    pvt[i]=i;
  for (i=0;i<n-1;i++) {
    maxind=i;
    max=fabs(A[i*n+i]);
    for (j=i+1;j<n;j++) {
      if (fabs(A[j*n+i])>max) {
        maxind = j;
        max = fabs(A[j*n+i]);
      }
    }
    if (max==0) {
      printf("LU factoration failed: coefficient matrix is singular\n\n");
      return -1;
    } else {
      if (maxind!=i) {
        temps=pvt[i];
        pvt[i]=pvt[maxind];
        pvt[maxind]=temps;
        memcpy(tempv,&A[i],sizeof(double)*n);
        memcpy(&A[i],&A[maxind],sizeof(double)*n);
        memcpy(&A[maxind],tempv,sizeof(double)*n); 
      }
    }
    for (j=i+1;j<n;j++) {
      A[j*n+i]=A[j*n+i]/A[i*n+i];
      for (k=i+1;k<n;k++)
        A[j*n+k]=A[j*n+k]-A[j*n+i]*A[i*n+k];
    }
  }
  printf("LU:\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      printf("%f ",A[i*n+j]);
    printf("\n");
  }
  y[0]=b[pvt[0]];
  for (i=1;i<n;i++) {
    sum=0;
    for (j=0;j<i;j++)
      sum+=y[j]*A[i*n+j];
    y[i]=b[pvt[i]]-sum;
  }
  x[n-1]=y[n-1]/A[(n-1)*n+n-1];
  for (i=n-1;i>-1;i--) {
    sum=0;
    for (j=i+1;j<n;j++)
      sum+=x[j]*A[i*n+j];
    x[i]=(y[i]-sum)/A[i*n+i];
  }
  printf("My result:\n");
  for (i=1;i<n;i++) 
    printf("%f ",x[i]);
  printf("\n");
  printf("Cost %.2f seconds by my approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  lapack_int *ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  double temp;
  
  start=clock();
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A_bak, n, ipiv);
  printf("LU:\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      printf("%f ",A[i*n+j]);
    printf("\n");
  }
  for (i=n-1;i>-1;i--) {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      temp=b[i];
      b[i]=b[ipiv[i]];
      b[ipiv[i]]=temp;
    }
  }
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, 1, 1, A_bak, n, b, 1);
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, A_bak, n, b, 1);
  printf("LA result:\n");
  for (i=1;i<n;i++) 
    printf("%f ",b[i]);
  printf("\n");
  printf("Cost %.2f seconds by LAPACKE's approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  if (test) {
    double max_diff=0, cur_diff;
    for (i=0;i<n;i++) {
      cur_diff=fabs(x[i]-b[i]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
    printf("The maximum difference between LAPACKE's approach and mine is %f.\n", max_diff);
  }
  
  printf("\n");

  free(A);
  free(A_bak);
  free(b);
  free(tempv);
  free(x);
  free(y);
  free(pvt);
  free(ipiv);

  return 0;
}
