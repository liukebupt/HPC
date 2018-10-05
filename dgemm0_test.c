#include <time.h>                                                                                                        
#include <stdio.h>                                                                                                       
#include <stdlib.h>                                                                                                      
#define drand() (double)rand()/RAND_MAX //return a random double number between 0 and 1.  
#define n 1024
                                                                                                                         
int main(int argc, char* argv[]) {       

  printf("Testing dgemm0 with n=%d.\n",n);                                                                               
                                                                                                                         
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
                                                                                                                         
  for (i=0;i<n;i++)                                                                                                      
    for (j=0;j<n;j++)                                                                                                    
      for (k=0;k<n;k++)                                                                                                  
        c[i*n+j] += a[i*n+k] * b[k*n+j];                                                                                 
                                                                                                                         
  printf("time passed %.2f seconds.\n",(double)(clock()-start)/CLOCKS_PER_SEC);                                          
                                                                                                                         
  free(a);                                                                                                               
  free(b);                                                                                                               
  free(c);                                                                                                               
                                                                                                                         
  return 0;                                                                                                              
}  
