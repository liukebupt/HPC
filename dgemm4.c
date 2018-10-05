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

  for (i=0;i<n;i+=2)
    for (j=0;j<n;j+=2) {
      register double c00=0, c01=0, c10=0, c11=0;
      for (k=0;k<n;k++) {
        c00 += a[i*n+k]*b[k*n+j];
        c10 += a[(i+1)*n+k]*b[k*n+j];
        c01 += a[i*n+k]*b[k*n+j+1];
        c11 += a[(i+1)*n+k]*b[k*n+j+1];
      }
      c[i*n+j]=c00;
      c[i*n+j+1]=c01;
      c[(i+1)*n+j]=c10;
      c[(i+1)*n+j+1]=c11;
    }
    
  printf("time passed %.2f seconds.\n",(double)(clock()-start)/CLOCKS_PER_SEC);

  free(a);
  free(b);
  free(c);

  return 0;
}
