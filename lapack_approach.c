#include <stdio.h>
#include <lapacke.h>
#define drand() (double)rand()/RAND_MAX+1 //return a random double number between 1 and 2.

int main (int argc, const char * argv[]) {

  if (argc!=2) {
    printf("Invalid input!\n");
    return 0;
  }
  
  int n=atoi(argv[1]);
  
  double *A=(double *)malloc(sizeof(double)*n*n);
  double *B=(double *)malloc(sizeof(double)*n);
  
  int i, j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      A[i][j]=drand();
    }
    B[i]=drand();
  }  

  free(A);
  free(B);

  return 0;
}
