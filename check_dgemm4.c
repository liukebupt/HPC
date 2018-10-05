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

  for (i=0;i<n;i+=2)
    for (j=0;j<n;j+=2) {
      register double c00=0, c01=0, c10=0, c11=0;
      for (k=0;k<n;k++) {
        c00 += a[i*n+k]*b[k*n+j];
        c10 += a[(i+1)*n+k]*b[k*n+j];
        c01 += a[i*n+k]*b[k*n+j+1];
        c11 += a[(i+1)*n+k]*b[k*n+j+1];
      }
      c3[i*n+j]=c00;
      c3[i*n+j+1]=c10;
      c3[(i+1)*n+j]=c01;
      c3[(i+1)*n+j+1]=c11;
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
