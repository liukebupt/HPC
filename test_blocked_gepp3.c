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
  double *b_bak=(double *)malloc(sizeof(double)*n);
  int i, j, k, l, j1, k1;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
    b[i]=drand(); 
  }
  memcpy(A_bak,A,sizeof(double)*n*n);
  memcpy(b_bak,b,sizeof(double)*n);
  
  int temps, maxind, end;
  double max, sum;
  double *tempv = (double *)malloc(sizeof(double)*n);
  double *y = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*n);
  
  int *pvt = (int *)malloc(sizeof(int)*n);

  clock_t start=clock();
  ///*
  printf("Blocked input:\n\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A[i*n+j]);
    }
    printf("\n\n");
  }
  //*/
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
    for (j=end;j<n;j+=B)
      for (k=end;k<n;k+=B)
        for (j1=j;j1<j+B;j1++)
          for (k1=k;k1<k+B;k1++) {
            register double a=A[j1*n+k1];
            for (l=i;l<end;l++)
              a-=A[j1*n+l]*A[l*n+k1];
            A[j1*n+k1]=a;
          }
  }
  ///*
  printf("Blocked LU:\n\n");
   for (i=0;i<n;i++) {
     for (j=0;j<n;j++) {
       printf("%f\t", A[i*n+j]);
     }
     printf("%d\n\n", pvt[i]);
  }
  //*/
  y[0]=b[pvt[0]];
  for (i=1;i<n;i++) {
    sum=0;
    for (j=0;j<i;j++)
      sum+=y[j]*A[i*n+j];
    y[i]=b[pvt[i]]-sum;
  }
  /*
  printf("blocked y:\n");
  for (k=0;k<n;k++) {
    printf("%f\t", y[k]);
  }
  printf("\n");
  */
  x[n-1]=y[n-1]/A[(n-1)*n+n-1];
  for (i=n-1;i>-1;i--) {
    sum=0;
    for (j=i+1;j<n;j++)
      sum+=x[j]*A[i*n+j];
    x[i]=(y[i]-sum)/A[i*n+i];
  }
  /*
  printf("blocked x:\n");
  for (k=0;k<n;k++) {
    printf("%f\t", x[k]);
  }
  printf("\n");
  */
  lapack_int *ipiv = (lapack_int *)malloc(n*sizeof(lapack_int));
  double temp;
  
  start=clock();
  memcpy(A,A_bak,sizeof(double)*n*n);
  printf("LAPACK input:\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A[i*n+j]);
    }
    printf("%f\n", b_bak[i]);
  }
  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
  printf("LAPACK LU:\n\n");
   for (i=0;i<n;i++) {
     for (j=0;j<n;j++) {
       printf("%f\t", A[i*n+j]);
     }
     printf("%d\n\n", ipiv[i]);
  }
  for (i=n-1;i>-1;i--) {
    ipiv[i]--;
    if (ipiv[i]!=i) {
      temp=b_bak[i];
      b_bak[i]=b_bak[ipiv[i]];
      b_bak[ipiv[i]]=temp;
    }
  }
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, 1, 1, A, n, b_bak, 1);
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, 1, 1, A, n, b_bak, 1);
  //printf("Cost %.2f seconds by LAPACKE's approach.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
  
  if (test) {
    double max_diff=0, cur_diff;
    for (i=0;i<n;i++) {
      cur_diff=fabs(x[i]-b_bak[i]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
    printf("The maximum difference between LAPACKE's approach and mine is %.16f.\n", max_diff);
  }
  
  memcpy(A,A_bak,sizeof(double)*n*n);
  /*
  printf("Simple LU input:\n\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      printf("%f\t", A[i*n+j]);
    }
    printf("\n\n");
  }
  */
  
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
  /*
  printf("Simple LU:\n\n");
   for (i=0;i<n;i++) {
     for (j=0;j<n;j++) {
       printf("%f\t", A[i*n+j]);
     }
     printf("%d\n\n", pvt[i]);
  }
  */
  y[0]=b[pvt[0]];
  for (i=1;i<n;i++) {
    sum=0;
    for (j=0;j<i;j++)
      sum+=y[j]*A[i*n+j];
    y[i]=b[pvt[i]]-sum;
  }
  /*
  printf("simple y:\n");
  for (k=0;k<n;k++) {
    printf("%f\t", y[k]);
  }
  printf("\n");
  */
  x[n-1]=y[n-1]/A[(n-1)*n+n-1];
  for (i=n-1;i>-1;i--) {
    sum=0;
    for (j=i+1;j<n;j++)
      sum+=x[j]*A[i*n+j];
    x[i]=(y[i]-sum)/A[i*n+i];
  }
  /*
  printf("simple x:\n");
  for (k=0;k<n;k++) {
    printf("%f\t", x[k]);
  }
  printf("\n");
  */
  if (test) {
    double max_diff=0, cur_diff;
    for (i=0;i<n;i++) {
      cur_diff=fabs(x[i]-b_bak[i]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
    printf("The maximum difference between LAPACKE's approach and mine is %.16f.\n", max_diff);
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
