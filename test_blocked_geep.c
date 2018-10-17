#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define drand() (double)rand()/RAND_MAX*(-2)+1 //return a random double number between -1 and 1.

int main (int argc, const char * argv[]) {

  int n=10, B=2;
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  int i, j, k,l;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
  printf("input:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%f ",A[i*n+j]);
     printf("\n");
  }
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind;
  double max;
  double *tempv = (double *)malloc(sizeof(double)*n);
  
  int *pvt = (int *)malloc(sizeof(int)*n);
  for (i=0;i<n;i++)
    pvt[i]=i;
  for (i=0;i<2-1;i++) {
    maxind=i;
    max=fabs(A_bak[i*n+i]);
    for (j=i+1;j<n;j++) {
      if (fabs(A_bak[j*n+i])>max) {
        maxind = j;
        max = fabs(A_bak[j*n+i]);
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
        memcpy(tempv,&A_bak[i*n],sizeof(double)*n);
        memcpy(&A_bak[i*n],&A_bak[maxind],sizeof(double)*n);
        memcpy(&A_bak[maxind*n],tempv,sizeof(double)*n); 
      }
    }
    for (j=i+1;j<n;j++) {
      A_bak[j*n+i]=A_bak[j*n+i]/A_bak[i*n+i];
      for (k=i+1;k<n;k++)
        A_bak[j*n+k]=A_bak[j*n+k]-A_bak[j*n+i]*A_bak[i*n+k];
    }
  }
  printf("Simple LU:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%f ",A_bak[i*n+j]);
     printf("  %d\n", pvt[i]);
  }
  
  int end;
  for (i=0;i<n;i++)
    pvt[i]=i;
  for (i=0;i<2;i+=B) {
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
          memcpy(&A[j*n],&A[maxind],sizeof(double)*n);
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
  
  printf("Blocked LU:\n");
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
