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
  int i, j, k, l, j1, k1;
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
    for (j=end;j<n;j+=B)
      for (k=end;k<n;k+=B)
        for (j1=j;j1<j+B;j1+=8)
          for (k1=k;k1<k+B;k1+=8) {
            register double c00=A[j1*n+k1], c01=A[j1*n+k1+1], c02=A[j1*n+k1+2], c03=A[j1*n+k1+3];
            register double c04=A[j1*n+k1+4], c05=A[j1*n+k1+5], c06=A[j1*n+k1+6], c07=A[j1*n+k1+7];
            register double c10=A[j1*n+k1+n], c11=A[j1*n+k1+1+n], c12=A[j1*n+k1+2+n], c13=A[j1*n+k1+3+n];
            register double c14=A[j1*n+k1+n+4], c15=A[j1*n+k1+n+5], c16=A[j1*n+k1+n+6], c17=A[j1*n+k1+n+7];
            register double c20=A[j1*n+k1+2*n], c21=A[j1*n+k1+1+2*n], c22=A[j1*n+k1+2+2*n], c23=A[j1*n+k1+3+2*n];
            register double c24=A[j1*n+k1+2*n+4], c25=A[j1*n+k1+2*n+5], c26=A[j1*n+k1+2*n+6], c27=A[j1*n+k1+2*n+7];
            register double c30=A[j1*n+k1+3*n], c31=A[j1*n+k1+1+3*n], c32=A[j1*n+k1+2+3*n], c33=A[j1*n+k1+3+3*n];
            register double c34=A[j1*n+k1+3*n+4], c35=A[j1*n+k1+3*n+5], c36=A[j1*n+k1+3*n+6], c37=A[j1*n+k1+3*n+7];
            register double c40=A[j1*n+k1+4*n], c41=A[j1*n+k1+1+4*n], c42=A[j1*n+k1+2+4*n], c43=A[j1*n+k1+3+4*n];
            register double c44=A[j1*n+k1+4*n+4], c45=A[j1*n+k1+4*n+5], c46=A[j1*n+k1+4*n+6], c47=A[j1*n+k1+4*n+7];
            register double c50=A[j1*n+k1+5*n], c51=A[j1*n+k1+1+5*n], c52=A[j1*n+k1+2+5*n], c53=A[j1*n+k1+3+5*n];
            register double c54=A[j1*n+k1+5*n+4], c55=A[j1*n+k1+5*n+5], c56=A[j1*n+k1+5*n+6], c57=A[j1*n+k1+5*n+7];
            register double c60=A[j1*n+k1+6*n], c61=A[j1*n+k1+1+6*n], c62=A[j1*n+k1+2+6*n], c63=A[j1*n+k1+3+6*n];
            register double c64=A[j1*n+k1+6*n+4], c65=A[j1*n+k1+6*n+5], c66=A[j1*n+k1+6*n+6], c67=A[j1*n+k1+6*n+7];
            register double c70=A[j1*n+k1+7*n], c71=A[j1*n+k1+1+7*n], c72=A[j1*n+k1+2+7*n], c73=A[j1*n+k1+3+7*n];
            register double c74=A[j1*n+k1+7*n+4], c75=A[j1*n+k1+7*n+5], c76=A[j1*n+k1+7*n+6], c77=A[j1*n+k1+7*n+7];
            for (l=i;l<end;l++) {
              register double a0=A[j1*n+l], a1=A[j1*n+l+n], a2=A[j1*n+l+2*n], a3=A[j1*n+l+3*n];
              register double a4=A[j1*n+l+4*n], a5=A[j1*n+l+5*n], a6=A[j1*n+l+6*n], a7=A[j1*n+l+7*n];
              register double b0=A[l*n+k1], b1=A[l*n+k1+1], b2=A[l*n+k1+2], b3=A[l*n+k1+3];
              register double b4=A[l*n+k1+4], b5=A[l*n+k1+5], b6=A[l*n+k1+6], b7=A[l*n+k1+7];
              c00-=a0*b0;
              c01-=a0*b1;
              c02-=a0*b2;
              c03-=a0*b3;
              c04-=a0*b4;
              c05-=a0*b5;
              c06-=a0*b6;
              c07-=a0*b7;
              c10-=a1*b0;
              c11-=a1*b1;
              c12-=a1*b2;
              c13-=a1*b3;
              c14-=a1*b4;
              c15-=a1*b5;
              c16-=a1*b6;
              c17-=a1*b7;
              c20-=a2*b0;
              c21-=a2*b1;
              c22-=a2*b2;
              c23-=a2*b3;
              c24-=a2*b4;
              c25-=a2*b5;
              c26-=a2*b6;
              c27-=a2*b7;
              c30-=a3*b0;
              c31-=a3*b1;
              c32-=a3*b2;
              c33-=a3*b3;
              c34-=a3*b4;
              c35-=a3*b5;
              c36-=a3*b6;
              c37-=a3*b7;
              c40-=a4*b0;
              c41-=a4*b1;
              c42-=a4*b2;
              c43-=a4*b3;
              c44-=a4*b4;
              c45-=a4*b5;
              c46-=a4*b6;
              c47-=a4*b7;
              c50-=a5*b0;
              c51-=a5*b1;
              c52-=a5*b2;
              c53-=a5*b3;
              c54-=a5*b4;
              c55-=a5*b5;
              c56-=a5*b6;
              c57-=a5*b7;
              c60-=a6*b0;
              c61-=a6*b1;
              c62-=a6*b2;
              c63-=a6*b3;
              c64-=a6*b4;
              c65-=a6*b5;
              c66-=a6*b6;
              c67-=a6*b7;
              c70-=a7*b0;
              c71-=a7*b1;
              c72-=a7*b2;
              c73-=a7*b3;
              c74-=a7*b4;
              c75-=a7*b5;
              c76-=a7*b6;
              c77-=a7*b7;
            }
            A[j1*n+k1]=c00;
            A[j1*n+k1+1]=c01;
            A[j1*n+k1+2]=c02;
            A[j1*n+k1+3]=c03;
            A[j1*n+k1+4]=c04;
            A[j1*n+k1+5]=c05;
            A[j1*n+k1+6]=c06;
            A[j1*n+k1+7]=c07;
            A[j1*n+k1+n]=c10;
            A[j1*n+k1+n+1]=c11;
            A[j1*n+k1+n+2]=c12;
            A[j1*n+k1+n+3]=c13;
            A[j1*n+k1+n+4]=c14;
            A[j1*n+k1+n+5]=c15;
            A[j1*n+k1+n+6]=c16;
            A[j1*n+k1+n+7]=c17;
            A[j1*n+k1+2*n]=c20;
            A[j1*n+k1+2*n+1]=c21;
            A[j1*n+k1+2*n+2]=c22;
            A[j1*n+k1+2*n+3]=c23;
            A[j1*n+k1+2*n+4]=c24;
            A[j1*n+k1+2*n+5]=c25;
            A[j1*n+k1+2*n+6]=c26;
            A[j1*n+k1+2*n+7]=c27;
            A[j1*n+k1+3*n]=c30;
            A[j1*n+k1+3*n+1]=c31;
            A[j1*n+k1+3*n+2]=c32;
            A[j1*n+k1+3*n+3]=c33;
            A[j1*n+k1+3*n+4]=c34;
            A[j1*n+k1+3*n+5]=c35;
            A[j1*n+k1+3*n+6]=c36;
            A[j1*n+k1+3*n+7]=c37;
            A[j1*n+k1+4*n]=c40;
            A[j1*n+k1+4*n+1]=c41;
            A[j1*n+k1+4*n+2]=c42;
            A[j1*n+k1+4*n+3]=c43;
            A[j1*n+k1+4*n+4]=c44;
            A[j1*n+k1+4*n+5]=c45;
            A[j1*n+k1+4*n+6]=c46;
            A[j1*n+k1+4*n+7]=c47;
            A[j1*n+k1+5*n]=c50;
            A[j1*n+k1+5*n+1]=c51;
            A[j1*n+k1+5*n+2]=c52;
            A[j1*n+k1+5*n+3]=c53;
            A[j1*n+k1+5*n+4]=c54;
            A[j1*n+k1+5*n+5]=c55;
            A[j1*n+k1+5*n+6]=c56;
            A[j1*n+k1+5*n+7]=c57;
            A[j1*n+k1+6*n]=c60;
            A[j1*n+k1+6*n+1]=c61;
            A[j1*n+k1+6*n+2]=c62;
            A[j1*n+k1+6*n+3]=c63;
            A[j1*n+k1+6*n+4]=c64;
            A[j1*n+k1+6*n+5]=c65;
            A[j1*n+k1+6*n+6]=c66;
            A[j1*n+k1+6*n+7]=c67;
            A[j1*n+k1+7*n]=c70;
            A[j1*n+k1+7*n+1]=c71;
            A[j1*n+k1+7*n+2]=c72;
            A[j1*n+k1+7*n+3]=c73;
            A[j1*n+k1+7*n+4]=c74;
            A[j1*n+k1+7*n+5]=c75;
            A[j1*n+k1+7*n+6]=c76;
            A[j1*n+k1+7*n+7]=c77;
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
