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
  register int n=atoi(argv[1]);
  register int B=atoi(argv[2]);
  printf("Testing blocked gepp with n=%d, B=%d.\n", n, B);
  
  bool test=false;
  if (argv[3][0]=='T')
    test=true;
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n);
  double *A_bak=(double *)malloc(sizeof(double)*n*n);
  register int i, j, k, l, j1, k1;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      A[i*n+j]=drand();
    b[i]=drand(); 
  }
  memcpy(A_bak,A,sizeof(double)*n*n);
  
  register int temps, maxind, end;
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
      register int p=j*n+j;
      max=fabs(A[p]);
      for (k=j+1;k<n;k++) {
        p+=n;
        if (fabs(A[p])>max) {
          maxind = k;
          max = fabs(A[p]);
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
      p=j*n+j;
      register int p1=p;       //k*n+j
      register double a1=A[p1];
      for (k=j+1;k<n;k++) {
        p1+=n;
        A[p1]/=a1;
        register double a=A[p1];
        register int p2=p1;        //k*n+l
        register int p3=p;     //j*n+l
        for (l=j+1;l<end;l++) {
          p3+=1;
          p2+=1;
          A[p2]-=a*A[p3];
        }
      }
    }
    
    register int p4=i*n
    register int p1=p4+i;         //k*n+i
    for (l=end;l<n;l++) {
      register int p=p4+l;         //k*n+l
      for (k=i+1;k<end;k++) {
        p1+=n;
        p+=n;
        register double a=A[p];
        register int p2=p1;     //k*n+j
        register int p3=p4+l;     //j*n+l
        for (j=i;j<k;j++) {
          a-=A[p2]*A[p3];
          p2+=1;
          p3+=n;
        }
        A[p]=a;
      }
      
    }
    
    for (j=end;j<n;j+=B)
      for (k=end;k<n;k+=B)
        for (j1=j;j1<j+B;j1+=4)
          for (k1=k;k1<k+B;k1+=4) {
            register double c00=A[j1*n+k1], c01=A[j1*n+k1+1], c02=A[j1*n+k1+2], c03=A[j1*n+k1+3];
            register double c10=A[j1*n+k1+n], c11=A[j1*n+k1+1+n], c12=A[j1*n+k1+2+n], c13=A[j1*n+k1+3+n];
            register double a0=A[j1*n+l], a1=A[j1*n+l+n], b0=A[l*n+k1], b1=A[l*n+k1+1];
            register double b2=A[l*n+k1+2], b3=A[l*n+k1+3];
            for (l=i;l<end;l++) {
              c00-=a0*b0;
              c01-=a0*b1;
              c02-=a0*b2;
              c03-=a0*b3;
              c10-=a1*b0;
              c11-=a1*b1;
              c12-=a1*b2;
              c13-=a1*b3;
            }
            A[j1*n+k1]=c00;
            A[j1*n+k1+1]=c01;
            A[j1*n+k1+2]=c02;
            A[j1*n+k1+3]=c03;
            A[j1*n+k1+n]=c10;
            A[j1*n+k1+n+1]=c11;
            A[j1*n+k1+n+2]=c12;
            A[j1*n+k1+n+3]=c13;
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
