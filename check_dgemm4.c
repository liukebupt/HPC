#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#define drand() (double)rand()/RAND_MAX //return a random double number between 0 and 1.

int main(int argc, char* argv[]) {
 
  if (argc!=2) {
    printf("Invalid input!\n");
    return 0;
  }

  //printf("RAND_MAX=%d, CLOCKS_PER_SEC=%d.\n",RAND_MAX,CLOCKS_PER_SEC);

  int n=atoi(argv[1]);

  printf("Checking dgemm4 with n=%d.\n",n);

  double *a=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n*n);
  double *c0=(double *)calloc(sizeof(double),n*n);
  double *c3=(double *)calloc(sizeof(double),n*n);

  int i;

  for (i=0; i<n*n; i++) {
    a[i]=drand();
    b[i]=drand();
  }
  
  int j,k;

  for (i=0;i<n;i+=2) {
    register int pc1=i*n;
    for (j=0;j<n;j+=4) {
      register int pc2=pc1+n, pa1=i*n, pb=j;
      register double c00=0, c01=0, c02=0, c03=0, c10=0, c11=0, c12=0, c13=0;
      for (k=0;k<n;k++) {
	       register int pa2=pa1+n;
	       register double a0=a[pa1], a1=a[pa2], b0=b[pb], b1=b[pb+1], b2=b[pb+2], b3=b[pb+3];
        c00 += a0*b0;
        c01 += a0*b1;
        c02 += a0*b2;
        c03 += a0*b3;
        c10 += a1*b0;
        c11 += a1*b1;
        c12 += a1*b2;
        c13 += a1*b3;
	       pa1+=1;
	       pb+=n;
      }
      c3[pc1]=c00;
      c3[pc1+1]=c01;
      c3[pc1+2]=c02;
      c3[pc1+3]=c03;
      c3[pc2]=c10;
      c3[pc2+1]=c11;
      c3[pc2+2]=c12;
      c3[pc2+3]=c13;
      pc1+=1;
    }
  }
  
  double max_diff=0,cur_diff;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++) {
      for (k=0;k<n;k++)
        c0[i*n+j]+=a[i*n+k]*b[k*n+j];
      cur_diff=abs(c0[i*n+j]-c3[i*n+j]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }

  printf("The maximum difference between dgemm0 and dgemm4 is %f.\n", max_diff);

  free(a);
  free(b);
  free(c0);
  free(c3);

  return 0;
}
