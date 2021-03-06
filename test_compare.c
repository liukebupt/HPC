#include <stdio.h>
#include <lapacke.h>
#include <time.h>
#include <stdbool.h>
#include <cblas.h>
#include <math.h>
#include <string.h>

#define drand() (double)rand()/RAND_MAX*(-2)+1 //return a random double number between -1 and 1.

int main (int argc, const char * argv[]) {

  if (argc!=3) {
    printf("Invalid input!\n");
    return 0;
  }
  int n=atoi(argv[1]);
  printf("n=%d. \n",n);
 
  bool test=false;
  if (argv[2][0]=='T')
    test=true;

  double *A=(double *)malloc(sizeof(double)*n*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  int i, j, k;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
    b[i]=drand();
  }  
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind;
  double sum, max;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  
  printf("input:\n\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A[i*n+j]);
    }
    printf("%f\n\n", b[i]);
  }
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
        memcpy(tempv,&A[i*n],sizeof(double)*n);
        memcpy(&A[i*n],&A[maxind*n],sizeof(double)*n);
        memcpy(&A[maxind*n],tempv,sizeof(double)*n); 
      }
    }
    for (j=i+1;j<n;j++) {
      A[j*n+i]=A[j*n+i]/A[i*n+i];
      for (k=i+1;k<n;k++)
        A[j*n+k]=A[j*n+k]-A[j*n+i]*A[i*n+k];
    }
  }
  printf("LU:\n\n");
  for (k=0;k<n;k++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A[k*n+j]);
    }
    printf("%d\n\n", pvt[k]);
  }
  y[0]=b[pvt[0]];
  for (i=1;i<n;i++) {
    sum=0;
    for (j=0;j<i;j++)
      sum+=y[j]*A[i*n+j];
    y[i]=b[pvt[i]]-sum;
  }
  printf("y:\n\n");
  for (k=0;k<n;k++) {
    printf("%f\t", y[k]);
  }
  printf("\n\n");
  x[n-1]=y[n-1]/A[(n-1)*n+n-1];
  for (i=n-1;i>-1;i--) {
    sum=0;
    for (j=i+1;j<n;j++)
      sum+=x[j]*A[i*n+j];
    x[i]=(y[i]-sum)/A[i*n+i];
  }
  printf("x:\n\n");
  for (k=0;k<n;k++) {
    printf("%f\t", x[k]);
  }
  printf("Cost %.2f seconds by my approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  lapack_int *ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  double temp;
  
  start=clock();
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A_bak, n, ipiv);
  printf("LAPACK LU:\n\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A_bak[i*n+j]);
    }
    printf("%d\n\n", pvt[i]);
  }
  for (i=0;i<n;i--) {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      temp=b[i];
      b[i]=b[ipiv[i]];
      b[ipiv[i]]=temp;
    }
  }
  printf("LAPACK b:\n\n");
  for (k=0;k<n;k++) {
    printf("%f\t", b[k]);
  }
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, 1, 1, A_bak, n, b, 1);
  printf("LAPACK y:\n\n");
  for (k=0;k<n;k++) {
    printf("%f\t", b[k]);
  }
  printf("\n\n");
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, A_bak, n, b, 1);
  printf("LAPACK x:\n\n");
  for (k=0;k<n;k++) {
    printf("%f\t", b[k]);
  }
  printf("Cost %.2f seconds by LAPACKE's approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  if (test) {
    double max_diff=0, cur_diff;
    for (i=0;i<n;i++) {
      cur_diff=fabs(x[i]-b[i]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
    printf("The maximum difference between LAPACKE's approach and mine is %.16f.\n", max_diff);
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
