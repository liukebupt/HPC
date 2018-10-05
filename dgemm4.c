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

  register int n=atoi(argv[1]);

  printf("Testing dgemm3 with n=%d.\n",n);

  double *a=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n*n);
  double *c=(double *)calloc(sizeof(double),n*n);

  register int i;

  for (i=0; i<n*n; i++) {
    a[i]=drand();
    b[i]=drand();
  }
  
  clock_t start=clock();
  
  register int j,k;

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
      c[pc1]=c00;
      c[pc1+1]=c01;
      c[pc1+2]=c02;
      c[pc1+3]=c03;
      c[pc2]=c10;
      c[pc2+1]=c11;
      c[pc2+2]=c12;
      c[pc2+3]=c13;
      pc1+=4;
    }
  }

  printf("time passed %.2f seconds.\n",(double)(clock()-start)/CLOCKS_PER_SEC);

  free(a);
  free(b);
  free(c);

  return 0;
}
