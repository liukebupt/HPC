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

  int B=atoi(argv[1]);

  printf("Testing blocked version with B=%d.\n",B);

  double *a=(double *)malloc(sizeof(double)*n*n);
  double *b=(double *)malloc(sizeof(double)*n*n);
  double *c_backup=(double *)malloc(sizeof(double)*n*n);

  int i;

  for (i=0; i<n*n; i++) {
    a[i]=drand();
    b[i]=drand();
    c_backup[i]=drand();
  }
  
  double *c, *c_answer；
  memcpy(c_answer,c_backup,sizeof(double)*n*n);
  
  int j, k;
  for (i=0;i<n;i++)
    for (j=0;j<n;j++) 
      for (k=0;k<n;k++)
        c_answer[i*n+j]+=a[i*n+k]*b[k*n+j];
        
  double max_diff, cur_diff;
  
  clock_t start;
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order ijk.\n");
  
  start=clock();
  for (i=0;i<n;i+=B)
    for (j=0;j<n;j+=B)
      for (k=0;k<n;k+=B)
        for (i1=i;i1<B;i1++)
          for (j1=j;j1<B;j1++) {
            register double r=c[i1*n+j1];
            for (k1=k;k1<B;k1++)
              r += a[i1*n+k1]*b[k1*n+j1];
            c[i1*n+j1]=r;
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
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order ikj.\n");
  
  start=clock();
  for (i=0;i<n;i+=B)
    for (k=0;k<n;k+=B)
      for (j=0;j<n;j+=B)
        for (i1=i;i1<B;i1++)
          for (k1=k;k1<B;k1++) {
            register double r=a[i1*n+k1];
            for (j1=j;j1<B;j1++)
              c[i1*n+j1] += r*b[k1*n+j1];
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
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order jik.\n");
  
  start=clock();
  for (j=0;j<n;j+=B)
    for (i=0;i<n;i+=B)
      for (k=0;k<n;k+=B)
        for (j1=j;j1<B;j1++)
          for (i1=i;i1<B;i1++) {
            register double r=c[i1*n+j1];
            for (k1=k;k1<B;k1++)
              r += a[i1*n+k1]*b[k1*n+j1];
            c[i1*n+j1]=r;
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
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order jki.\n");
  
  start=clock();
  for (j=0;j<n;j+=B)
    for (k=0;k<n;k+=B)
      for (i=0;i<n;i+=B)
        for (j1=j;j1<B;j1++)
          for (k1=k;k1<B;k1++) {
            register double r=b[k1*n+j1];
            for (i1=i;i1<B;i1++)
              c[i1*n+j1] += a[i1*n+k1]*r;
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
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order kij.\n");
  
  start=clock();
  for (k=0;k<n;k+=B)
    for (i=0;i<n;i+=B)
      for (j=0;j<n;j+=B)
        for (k1=k;k1<B;k1++)
          for (i1=i;i1<B;i1++) {
            register double r=a[i1*n+k1];
            for (j1=j;j1<B;j1++)
              c[i1*n+j1] += r*b[k1*n+j1];
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
  
  memcpy(c,c_backup,sizeof(double)*n*n);
  
  printf("Running blocked version for order kji.\n");
  
  start=clock();
  for (k=0;k<n;k+=B)
    for (j=0;j<n;j+=B)
      for (i=0;i<n;i+=B)
        for (k1=k;k1<B;k1++)
          for (j1=j;j1<B;j1++) {
            register double r=b[k1*n+j1];
            for (i1=i;i1<B;i1++)
              c[i1*n+j1] += a[i1*n+k1]*r;
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
  free(c_backup);
  free(c_answer)

  return 0;
}
