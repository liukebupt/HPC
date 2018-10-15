#include <stdio.h>
#include <lapacke.h>
#include <time.h>
#include <stdbool.h>
#include <cblas.h>
#include <string.h>

#define drand() (double)rand()/RAND_MAX+1 //return a random double number between 1 and 2.

int main (int argc, const char * argv[]) {

  if (argc!=3) {
    printf("Invalid input!\n");
    return 0;
  }
  int n=atoi(argv[1]);
  printf("n=.\n",n);
 
  bool test=false;
  if (argv[2]=="True")
    bool test=true;

  double *A=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  int i, j, k;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      A[i*n+j]=drand();
    }
    b[i]=drand();
  }  
  
  int temps, maxind, max;
  double sum;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  clock_t start=clock();
  int *pvt = (int *)malloc(sizeof(int)*n);
  for (i=0;i<n;i++)
    pvt[i]=i;
  for (i=0;i<n-1;i++) {
    maxind=i;
    max=abs(A[i*n+i]);
    for (j=i+1;j<n;j++) {
      if (abs(A[j*n+i])>max) {
        maxind = j;
        max = abs(A[j*n+i]);
      }
    }
    if (max==0) {
      printf("LUfactoration failed: coefficient matrix is singular");
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
  printf("Cost %.2f seconds by my approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  lapack_int *ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  double temp;
  
  start=clock();
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
  for (i=0;i<n;i++)
  {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      temp=b[i];
      b[i]=b[ipiv[i]];
      b[ipiv[i]]=temp;
    }
  }
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, 1, 1, A, n, b, 1);
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, A, n, b, 1);
  printf("Cost %.2f seconds by LAPACKE's approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  if (test) {
    double max_diff=0, cur_diff;
    for (i=0;i<n;i++) {
      cur_diff=abs(x[i]-b[i]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
    printf("The maximum difference between LAPACKE's approach and mine is %f.\n", max_diff);
  }
  
  printf("\n");

  free(A);
  free(b);
  free(tempv);
  free(x);
  free(y);
  free(pvt);
  free(ipiv);

  return 0;
}
