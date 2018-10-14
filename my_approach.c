#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, const char * argv[]) {

  double A[2][2] = {3,4,5,9}, b[] = {5,3}, *tempv;

  int i, j, k, n=2, *pvt, temps;
  
  pvt = (int *)malloc(sizeof(int)*n);
  tempv = (double *)malloc(sizeof(double)*n);
  
  for (i=0;i<n-1;i++) {
    pvt[i]=i;
  }
  
  for (i=0;i<n-1;i++) {
    int maxind=i, max=abs(A[i][i]);
    for (j=i+1;j<n;j++) {
      if (abs(A[j][i])>max) {
        maxind = j;
        max = abs(A[j][i]);
      }
    }
    if (max==0) {
      printf("LUfactoration failed: coefficient matrix is singular");
      return -1;
    } else {
      if (maxind!=i) {
        temps=pvt[i];
        pvt[i]=pvt[maxind];
        pvt[maxind]=temps;
        memcpy(tempv,A[i],sizeof(double)*n);
        memcpy(A[i],A[maxind],sizeof(double)*n);
        memcpy(A[maxind],tempv,sizeof(double)*n); 
        printf("After:\n");
        for(i=0;i<n;i++)
        {
           for(j=0;j<n;j++)
              printf("%lf ",A[i][j]);
           printf("  %6d \n", pvt[i]);
        }  
      }
    }
    for (j=i+1;j<n;j++) {
      A[j][i]=A[j][i]/A[i][i];
      for (k=i+1;k<n;k++)
        A[j][k]=A[j][k]-A[j][i]*A[i][k];
    }
  }
  
  printf("LU:\n");
  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
        printf("%lf ",A[i][j]);
     printf("  %6d \n", pvt[i]);
  }
  
  double y[2];
  
  y[0]=b[pvt[0]];
  for (i=1;i<n;i++) {
    double sum=0;
    for (j=0;j<i-1;j++)
      sum+=y[j]*A[i][j];
    y[i]=b[pvt[i]]-sum;
  }
  printf("y:\n");
  for(i=0;i<n;i++)
     printf("  %lf \n", y[i]);

  return 0;
}
