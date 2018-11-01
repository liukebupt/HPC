#include <mpi.h>     
#include <math.h>    
#include <stdio.h>   
#include <stdlib.h>  
#include "MyMPI.h"   
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])  
{      
  long long int n, low_value, Blow_value, Bsize, size, proc0_size, i, first, prime, gsize;
  int id, p, index, count, global_count, j, B, k, sqn, psize, *newprimes, nprime; 
  char *marked, *pend, *primes;     
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
  gsize=n/2-1;
  B = atoi(argv[2]);
  
  low_value = 3 + 2*BLOCK_LOW(id,p,gsize);     
  size = BLOCK_SIZE(id,p,gsize);  
  proc0_size = gsize/p;  
  sqn=(int) sqrt((double) n);
  if ((3 + 2*proc0_size) < sqn) {   
    if (!id) printf ("Too many processes\n");  
    MPI_Finalize();
    exit (1);      
  }   
  //printf("id = %d p = %d n = %lld size = %" PRId64 "\n", id, p, n, size);     

  psize=sqn/2;
  primes = (char *) calloc (psize, sizeof(char)); 
  index = 0;
  prime = 3; 
  do {
    first = (prime * prime - 3)/2;
    for (i = first; i < psize; i += prime) primes[i] = 1;     
    while (primes[++index]);  
    //printf("%d\n", prime);  
    prime = 2*index + 3;      
  } while (prime * prime <= sqn);
  newprimes = (int *) malloc (psize*sizeof(int)); 
  nprime=0;
  for (j = 0; j < psize; j++) {
    if (!primes[j]) {
      newprimes[nprime++]=2*j + 3;
      //if (!id && nprime<100) printf("%d\n", newprimes[nprime-1]);
   }
  }
  //printf("%d aaaaa %d\n", primes[14], psize);

  marked = (char *) calloc (size, sizeof(char));    
  if (marked == NULL) {    
    printf ("Cannot allocate enough memory\n");
    MPI_Finalize();
    exit (1);      
  }   
  //if (id==7) printf("aaaa %d\n", high_value);
  for (k=0;k<size;k+=B) {
    Blow_value = low_value+2*k;
    //printf("%d  %d\n", k+B, size);
    Bsize = MIN(k+B, size);
    //for (i = 0; i < Bsize; i ++) printf("%d\n", i*2+Blow_value);  
    for (j = 0; j < nprime; j++) {
      prime = newprimes[j]; 
      //if (k==0 && !id) printf("%d\n", prime);
      if (prime * prime > Blow_value)      
      first = (prime * prime - low_value)/2;  
      else {  
        if (!(Blow_value % prime)) first = k;    
        else first = prime - (Blow_value*(prime +1)/2 % prime)+k;
      }
      for (i = first; i < Bsize; i += prime) { marked[i] = 1;  
       // printf("%d\n", Blow_value+i*2);
      } 
      //if (!id) printf("first %d prime %d Bsize %d\n", first*2+low_value, prime, Bsize); 
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
