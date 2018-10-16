#include <stdio.h>
#include <lapacke.h>
//#include <time.h>
#include <stdbool.h>
#include <cblas.h>
#include <math.h>
#include <string.h>

#define drand() (double)rand()/RAND_MAX*(-2)+1 //return a random double number between -1 and 1.

int main (int argc, const char * argv[]) {

  int n=10;
  bool test=true
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  int i, j, k;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind;
  double sum, max;
  double *tempv = (double *)malloc(sizeof(double)*n);
  
  //clock_t start=clock();
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
        memcpy(&A[i*n],&A[maxind],sizeof(double)*n);
        memcpy(&A[maxind*n],tempv,sizeof(double)*n); 
      }
    }
    for (j=i+1;j<n;j++) {
      A[j*n+i]=A[j*n+i]/A[i*n+i];
      for (k=i+1;k<n;k++)
        A[j*n+k]=A[j*n+k]-A[j*n+i]*A[i*n+k];
    }
  }
  printf("Simple LU:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%f ",A[i*n+j]);
     printf("  %d\n", pvt[i]);
  }
  
  
  

  free(A);
  free(A_bak);
  free(tempv);
  free(pvt);

  return 0;
}
