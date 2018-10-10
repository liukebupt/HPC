#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#define drand() (double)rand()/RAND_MAX //return a random double number between 0 and 1.
#define n 2048

int main(int argc, char* argv[]) {
  if (argc!=2) {
    printf("Invalid input!\n");
    return 0;
  }

  register int B=atoi(argv[1]);

  printf("Testing blocked and register reused version with B=%d.\n",B);

  double *a=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n*n);
  double *c_backup=(double *)malloc(sizeof(double)*n*n);

  register int i;

  for (i=0; i<n*n; i++) {
    a[i]=drand();
    b[i]=drand();
    c[i]=drand();
  }
  
  double *c, *c_answer；
  memcpy(c_answer,c,sizeof(double)*n*n);
  
  register int j, k;
  for (i=0;i<n;i++)
    for (j=0;j<n;j++) 
      for (k=0;k<n;k++)
        c_answer[i*n+j]+=a[i*n+k]*b[k*n+j];
        
  double max_diff, cur_diff;
  
  clock_t start;
  
  register int i1,j1,k1
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order ijk.\n");
  
  start=clock();
  for (i=0;i<n;i+=B)
    for (j=0;j<n;j+=B)
      for (k=0;k<n;k+=B)
        for (i1=i;i1<B;i1+=2)
          for (j1=j;j1<B;j1+=4) {
            register int pc1=i1*n+j1, pc2=pc1+n;
            register double c00=c[pc1], c01=c[pc1+1], c02=c[pc1+2], c03=c[pc1+3], c10=c[pc2], c11=c[pc2+1], c12=c[pc2+2], c13=c[pc2+3];
            for (k1=k;k1<B;k1++) {
              register int pa1=i1*n+k1, pa2=pa1+n, pb=k1*n+j1;
              register double a0=a[pa1], a1=a[pa2], b0=b[pb], b1=b[pb+1], b2=b[pb+2], b3=b[pb+3];
              c00 += a0*b0;
              c01 += a0*b1;
              c02 += a0*b2;
              c03 += a0*b3;
              c10 += a1*b0;
              c11 += a1*b1;
              c12 += a1*b2;
              c13 += a1*b3;
            }
            c[pc1]=c00;
            c[pc1+1]=c01;
            c[pc1+2]=c02;
            c[pc1+3]=c03;
            c[pc2]=c10;
            c[pc2+1]=c11;
            c[pc2+2]=c12;
            c[pc2+3]=c13;
          }
      
  printf("time passed %.2f seconds.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
      
  max_diff=0;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++) {
      cur_diff=abs(c_answer[i*n+j]-c[i*n+j]);
      if (cur_diff>max_diff)
        max_diff=cur_diff;
    }
  
  printf("The maximum difference is %f.\n\n", max_diff);

  free(a);
  free(b);
  free(c);
  free(c_answer)

  return 0;
}
