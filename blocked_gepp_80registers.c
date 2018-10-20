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
  register double max, sum;
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
      register int p=j*n+j;  //k*n+j
      max=A[p];
      if (max<0)
        max*=(-1);
      for (k=j+1;k<n;k++) {
        p+=n;
        register double abs=A[p];  //abs(A[k*n+j])
        if (abs<0)
          abs*=(-1);
        if (abs>max) {
          maxind = k;
          max = abs;
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
    register int p4=i*n;
    for (l=end;l<n;l++) {
      register int p=p4+l;         //k*n+l
      register int p1=p4+i;         //k*n+i
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
        for (j1=j;j1<j+B;j1+=8) {
          register int p5=j1*n; //j1*n
          register int p1=p5+k;    //j1*n+k1
          for (k1=k;k1<k+B;k1+=8) {
            register double c00=A[p1], c01=A[p1+1], c02=A[p1+2], c03=A[p1+3];
            register double c04=A[p1+4], c05=A[p1+5], c06=A[p1+6], c07=A[p1+7];
            register int p2=p1+n; //j1*n+k1+n
            register double c10=A[p2], c11=A[p2+1], c12=A[p2+2], c13=A[p2+3];
            register double c14=A[p2+4], c15=A[p2+5], c16=A[p2+6], c17=A[p2+7];
            p2+=n;   //j1*n+k1+2*n
            register double c20=A[p2], c21=A[p2+1], c22=A[p2+2], c23=A[p2+3];
            register double c24=A[p2+4], c25=A[p2+5], c26=A[p2+6], c27=A[p2+7];
            p2+=n;   //j1*n+k1+3*n
            register double c30=A[p2], c31=A[p2+1], c32=A[p2+2], c33=A[p2+3];
            register double c34=A[p2+4], c35=A[p2+5], c36=A[p2+6], c37=A[p2+7];
            p2+=n;   //j1*n+k1+4*n
            register double c40=A[p2], c41=A[p2+1], c42=A[p2+2], c43=A[p2+3];
            register double c44=A[p2+4], c45=A[p2+5], c46=A[p2+6], c47=A[p2+7];
            p2+=n;   //j1*n+k1+5*n
            register double c50=A[p2], c51=A[p2+1], c52=A[p2+2], c53=A[p2+3];
            register double c54=A[p2+4], c55=A[p2+5], c56=A[p2+6], c57=A[p2+7];
            p2+=n;   //j1*n+k1+6*n
            register double c60=A[p2], c61=A[p2+1], c62=A[p2+2], c63=A[p2+3];
            register double c64=A[p2+4], c65=A[p2+5], c66=A[p2+6], c67=A[p2+7];
            p2+=n;   //j1*n+k1+7*n
            register double c70=A[p2], c71=A[p2+1], c72=A[p2+2], c73=A[p2+3];
            register double c74=A[p2+4], c75=A[p2+5], c76=A[p2+6], c77=A[p2+7];
            register int p3=p5+i; //j1*n+l
            for (l=i;l<end;l++) {
              register double a0=A[p3], a1=A[p3+n], a2=A[p3+2*n], a3=A[p3+3*n];
              register double a4=A[p3+4*n], a5=A[p3+5*n], a6=A[p3+6*n], a7=A[p3+7*n];
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
              p3+=1;
            }
            A[p1]=c00;
            A[p1+1]=c01;
            A[p1+2]=c02;
            A[p1+3]=c03;
            A[p1+4]=c04;
            A[p1+5]=c05;
            A[p1+6]=c06;
            A[p1+7]=c07;
            p2=p1+n;
            A[p2]=c10;
            A[p2+1]=c11;
            A[p2+2]=c12;
            A[p2+3]=c13;
            A[p2+4]=c14;
            A[p2+5]=c15;
            A[p2+6]=c16;
            A[p2+7]=c17;
            p2+=n;
            A[p2]=c20;
            A[p2+1]=c21;
            A[p2+2]=c22;
            A[p2+3]=c23;
            A[p2+4]=c24;
            A[p2+5]=c25;
            A[p2+6]=c26;
            A[p2+7]=c27;
            p2+=n;
            A[p2]=c30;
            A[p2+1]=c31;
            A[p2+2]=c32;
            A[p2+3]=c33;
            A[p2+4]=c34;
            A[p2+5]=c35;
            A[p2+6]=c36;
            A[p2+7]=c37;
            p2+=n;
            A[p2]=c40;
            A[p2+1]=c41;
            A[p2+2]=c42;
            A[p2+3]=c43;
            A[p2+4]=c44;
            A[p2+5]=c45;
            A[p2+6]=c46;
            A[p2+7]=c47;
            p2+=n;
            A[p2]=c50;
            A[p2+1]=c51;
            A[p2+2]=c52;
            A[p2+3]=c53;
            A[p2+4]=c54;
            A[p2+5]=c55;
            A[p2+6]=c56;
            A[p2+7]=c57;
            p2+=n;
            A[p2]=c60;
            A[p2+1]=c61;
            A[p2+2]=c62;
            A[p2+3]=c63;
            A[p2+4]=c64;
            A[p2+5]=c65;
            A[p2+6]=c66;
            A[p2+7]=c67;
            p2+=n;
            A[p2]=c70;
            A[p2+1]=c71;
            A[p2+2]=c72;
            A[p2+3]=c73;
            A[p2+4]=c74;
            A[p2+5]=c75;
            A[p2+6]=c76;
            A[p2+7]=c77;
            p1+=8;
          }
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
