#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <time.h>
#include <stdbool.h>

#define drand() (double)rand()/RAND_MAX*(-2)+1 //return a random double number between -1 and 1.

int main (int argc, const char * argv[]) {
  
  if (argc!=4) {
    printf("Invalid input!\n");
    return 0;
  }
  int n=atoi(argv[1]);
  int B=atoi(argv[2]);
  printf("Testing blocked gepp with n=%d, B=%d.\n", n, B);
  
  bool test=false;
  if (argv[3][0]=='T')
    test=true;
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  int i, j, k, l;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
    b[i]=drand(); 
  }
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind, end;
  double max, sum;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  int *pvt = (int *)malloc(sizeof(int)*n);

  clock_t start=clock();
  for (i=0;i<n;i++)
    pvt[i]=i;
  for (i=0;i<n;i+=B) {
    end=i+B;
    for (j=i;j<end;j++) {
      maxind=j;
      max=fabs(A[j*n+j]);
      for (k=j+1;k<n;k++) {
        if (fabs(A[k*n+j])>max) {
          maxind = k;
          max = fabs(A[k*n+j]);
        }
      }
      if (max==0) {
        printf("LU factoration failed: coefficient matrix is singular\n\n");
        return -1;
      } else {
        if (maxind!=j) {
          temps=pvt[j];
          pvt[j]=pvt[maxind];
          pvt[maxind]=temps;
          memcpy(tempv,&A[j*n],sizeof(double)*n);
          memcpy(&A[j*n],&A[maxind*n],sizeof(double)*n);
          memcpy(&A[maxind*n],tempv,sizeof(double)*n); 
        }
      }
      for (k=j+1;k<n;k++) {
        A[k*n+j]=A[k*n+j]/A[j*n+j];
        for (l=j+1;l<end;l++)
          A[k*n+l]=A[k*n+l]-A[k*n+j]*A[j*n+l];
      }
    }
    for (j=i;j<end;j++)
      for (k=j+1;k<end;k++)
        for (l=end;l<n;l++)
          A[k*n+l]-=A[k*n+j]*A[j*n+l];
    for (j=end;j<n;j++)
      for (k=end;k<n;k++)
        for (l=i;l<end;l++)
          A[j*n+k]-=A[j*n+l]*A[l*n+k];
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
  
  if (test) {
    double temp;
    int *ipiv = (int *)malloc(sizeof(int)*n);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A_bak, n, ipiv);
    for (i=0;i<n;i++)
    {
      ipiv[i]--;
      if (ipiv[i]!=i) {
        temp=b[i];
        b[i]=b[ipiv[i]];
        b[ipiv[i]]=temp;
      }
    }
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, 1, 1, A_bak, n, b, 1);
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, A_bak, n, b, 1);
  
    double cur_diff, max_diff=0;
    for(i=0;i<n;i++) {
      cur_diff=fabs((b[i]-x[i])/b[i]);
      if (cur_diff > max_diff)
        max_diff=cur_diff;
    }
    printf("max relative difference is %.16f.\n", max_diff);
    free(ipiv);
  }
  
  free(A_bak);
  free(A);
  free(b);
  free(tempv);
  free(pvt);
  free(y);
  free(x);

  return 0;
}
