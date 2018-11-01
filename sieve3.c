#include <mpi.h>     
#include <math.h>    
#include <stdio.h>   
#include <stdlib.h>  
#include "MyMPI.h"   
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])  
{      
  long long int n, low_value, high_value, Bhigh_value, size, proc0_size, i, first, prime;
  int id, p, index, count, global_count, j, B, k; 
  char *marked, *pend;     
  double elapsed_time;     
  MPI_Init (&argc, &argv); 
  MPI_Barrier(MPI_COMM_WORLD);    
  elapsed_time = -MPI_Wtime();    
  MPI_Comm_rank (MPI_COMM_WORLD, &id);   
  //printf("process %d starts\n", id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);    
  if (argc != 3) {  
    if (!id) printf ("Command line: %s <m>\n", argv[0]);     
    MPI_Finalize(); exit (1);    
  }   
  n = strtoll(argv[1], &pend, 10);  
  B = atoi(argv[2]);
  
  high_value = 3 + 2*BLOCK_HIGH(id,p,n/2-1);    
  size = BLOCK_SIZE(id,p,n/2-1);  
  proc0_size = (n/2-1)/p;  
  if ((3 + 2*proc0_size) < (int) sqrt((double) n)) {   
    if (!id) printf ("Too many processes\n");  
    MPI_Finalize();
    exit (1);      
  }   
  //printf("id = %d p = %d n = %lld size = %" PRId64 "\n", id, p, n, size);     

  int sqn=sqrt(n);
  int psize=sqn/2-1;
  char *primes = (char *) malloc (psize); 
  for (i = 0; i < psize; i++) primes[i] = 0;  
  index = 0;
  prime = 3; 
  do {
    first = (prime * prime - 3)/2;
    for (i = first; i < psize; i += prime) primes[i] = 1;     
    while (primes[++index]);  
    //printf("%d\n", prime);  
    prime = 2*index + 3;      
  } while ((long long int)prime * prime <= sqn);
  int *newprimes = (int *) malloc (psize*sizeof(int)); 
  int nprime=0;
  for (i = 0; i < psize; i++) {
    if (!primes[i]) {
      newprimes[nprime++]=2*i + 3;
      //if (!id && nprime<100) printf("%d\n", newprimes[nprime-1]);
   }
  }

  marked = (char *) malloc (size);    
  if (marked == NULL) {    
    printf ("Cannot allocate enough memory\n");
    MPI_Finalize();
    exit (1);      
  }   
  for (i = 0; i < size; i++) marked[i] = 0;  
  
  for (k=0;k<size;k+=B) {
    low_value = 3 + 2*(BLOCK_LOW(id,p,n/2-1)+k);
    Bhigh_value = MIN(low_value+2*(B-1), high_value);
    for (j = 0; j < nprime; j++) {
      prime = newprimes[j]; 
      if (prime * prime > low_value)      
      first = (prime * prime - low_value)/2;  
      else {  
        if (!(low_value % prime)) first = 0;    
        else first = prime - (low_value*(prime +1)/2 % prime);
      }
      for (i = first; i < Bhigh_value; i += prime) marked[i] = 1;    
    }
  }
  count = 0; 
  for (i = 0; i < size; i++)      
    if (!marked[i]) { count++;   
//printf("%d\n",(i*2+low_value));
    }
  MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);      
  elapsed_time += MPI_Wtime();    
  if (!id) { 
    printf ("%lld primes are less than or equal to %lld\n", global_count+1, n);   
    printf ("Total elapsed time: %10.6f for %d cores blocksize %d.\n", elapsed_time, p, B);
  }   
  MPI_Finalize ();  
  return 0;  
}  
