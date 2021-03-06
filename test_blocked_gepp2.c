#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>

#define drand() ((double)rand()/RAND_MAX*(-2)+1) //return a random double number between -1 and 1.

int main (int argc, const char * argv[]) {
  
  if (argc!=4) {
    printf("Invalid input!\n");
    return 0;
  }

  int size=atoi(argv[1]);
  int n=atoi(argv[2]);
  int B=atoi(argv[3]);

  printf("Testing blocked gepp with size=%d, n=%d, B=%d.\n", size, n, B);
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  int i, j, k,l;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      A[i*n+j]=drand()*size;
    b[i]=drand()*size; 
  }
  printf("input:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%f\t",A[i*n+j]);
     printf("\n\n");
  }
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  int temps, maxind;
  double max, sum, temp;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  int *pvt = (int *)malloc(sizeof(int)*n);
  int *ipiv = (int *)malloc(sizeof(int)*n);
  
  int end;
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
      printf("Blocked LU after %d pivot:\n", j+1);
      for (k=0;k<n;k++) {
        for (l=0;l<n;l++)
          printf("%f\t", A[k*n+l]);
        printf("%d\n\n", pvt[k]);
      }
      for (k=j+1;k<n;k++) {
        A[k*n+j]=A[k*n+j]/A[j*n+j];
        for (l=j+1;l<end;l++)
          A[k*n+l]=A[k*n+l]-A[k*n+j]*A[j*n+l];
      }
      printf("Blocked LU after %d:\n", j+1);
      for (k=0;k<n;k++) {
        for (l=0;l<n;l++)
          printf("%f\t", A[k*n+l]);
        printf("%d\n\n", pvt[k]);
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
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++)
      printf("%f\t", A[i*n+j]);
    printf("%d\n\n", pvt[i]);
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
  
  memcpy(A,A_bak,sizeof(double)*n*n);
  
  printf("input:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%f\t",A[i*n+j]);
     printf("\n\n");
  }
  
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
        memcpy(tempv,&A[i*n],sizeof(double)*n);
        memcpy(&A[i*n],&A[maxind*n],sizeof(double)*n);
        memcpy(&A[maxind*n],tempv,sizeof(double)*n); 
      }
    }
    printf("Simple LU after %d pivot:\n", i+1);
    for (k=0;k<n;k++) {
      for (l=0;l<n;l++)
        printf("%f\t", A[k*n+l]);
      printf("%d\n\n", pvt[k]);
    }
    for (j=i+1;j<n;j++) {
      A[j*n+i]=A[j*n+i]/A[i*n+i];
      for (k=i+1;k<n;k++)
        A[j*n+k]=A[j*n+k]-A[j*n+i]*A[i*n+k];
    }
    printf("Simple LU after %d:\n", i+1);
    for (k=0;k<n;k++) {
      for (l=0;l<n;l++)
        printf("%f\t", A[k*n+l]);
      printf("%d\n\n", pvt[k]);
    }
  }
  printf("Simple LU:\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++)
      printf("%f\t", A[i*n+j]);
    printf("%d\n\n", pvt[i]);
  }
  
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

  free(A);
  free(A_bak);
  free(tempv);
  free(pvt);
  free(y);
  free(x);
  free(ipiv);

  return 0;
}
