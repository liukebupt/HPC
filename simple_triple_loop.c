#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#define drand() (double)rand()/RAND_MAX //return a random double number between 0 and 1.
#define n 2048

int main(int argc, char* argv[]) {

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
  
  printf("Running simple triple loop for order ijk.\n");
  
  start=clock();
  for (i=0;i<n;i++)
    for (j=0;j<n;j++) {
      register double r=c[i*n+j];
        for (k=0;k<n;k++)
          r += a[i*n+k]*b[k*n+j];
        c[i*n+j]=r;
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
  
  printf("Running simple triple loop for order ikj.\n");
  
  start=clock();
  for (i=0;i<n;i++)
    for (k=0;k<n;k++) {
      register double r=a[i*n+k];
        for (j=0;j<n;j++)
          c[i*n+j] += r*b[k*n+j];
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
  
  printf("Running simple triple loop for order jik.\n");
  
  start=clock();
  for (j=0;j<n;j++)
    for (i=0;i<n;i++) {
      register double r=c[i*n+j];
        for (k=0;k<n;k++)
          r += a[i*n+k]*b[k*n+j];
        c[i*n+j]=r;
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
  
  printf("Running simple triple loop for order jki.\n");
  
  start=clock();
  for (j=0;j<n;j++)
    for (k=0;k<n;k++) {
      register double r=b[k*n+j];
        for (i=0;i<n;i++)
          c[i*n+j] += a[i*n+k]*r;
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
  
  printf("Running simple triple loop for order kij.\n");
  
  start=clock();
  for (k=0;k<n;k++)
    for (i=0;i<n;i++) {
      register double r=a[i*n+k];
        for (j=0;j<n;j++)
          c[i*n+j] += r*b[k*n+j];
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
  
  printf("Running simple triple loop for order kji.\n");
  
  start=clock();
  for (k=0;k<n;k++)
    for (j=0;j<n;j++) {
      register double r=b[k*n+j];
        for (i=0;i<n;i++)
          c[i*n+j] += a[i*n+k]*r;
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
