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

  printf("Testing dgemm3 with n=%d.\n",n);

  double *a=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n*n);
  double *c=(double *)calloc(sizeof(double),n*n);

  int i;

  for (i=0; i<n*n; i++) {
    a[i]=drand();
    b[i]=drand();
  }
  
  clock_t start=clock();
  
  int j,k;

  register int pc1=0, pc2=n;

  for (i=0;i<n;i+=2) {
    for (j=0;j<n;j+=4) {
      register double c00=0, c01=0, c02=0, c03=0, c10=0, c11=0, c12=0, c13=0;
      register int pa1=i*n, pa2=pa1+n, pb=j;
      for (k=0;k<n;k++) {
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
	pa2+=1;
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
      pc2+=4;
    }
    pc1+=n;
    pc2+=n;
  }  

  printf("time passed %.2f seconds.\n",(double)(clock()-start)/CLOCKS_PER_SEC);

  free(a);
  free(b);
  free(c);

  return 0;
}
