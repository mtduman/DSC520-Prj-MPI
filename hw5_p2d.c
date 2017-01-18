#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h> // include mpi library header
#include <omp.h>
#include "secret_function.h"

void sample_rand(const double a, const double b ,const int dim, double *x) {
  for(int i=0;i<dim;++i) {
    double tmp = ((double) rand())/((double) RAND_MAX);
    tmp = (b-a)*tmp + a;
    x[i] = tmp;
  }
}

int main(int argc, char **argv)
{
  double start_time = omp_get_wtime();
// printf( "f( 512, 404.2319 ) = %f \n" ,secret_function( 512.0 , 404.2319) ); 

  //    const unsigned int N = 1.e7; // for testing 
  unsigned int N = atoi( argv[1] );
    
  int rank, size;

  MPI_Init (&argc, &argv);      // initializes MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); // get current MPI-process ID. O, 1, ...
  MPI_Comm_size (MPI_COMM_WORLD, &size); // get the total number of processes

  const int dim = 10;
  double x[dim];  // array of random numbers
  double min = 512., min_global;

  srand(time(NULL) * (int)rank); // each MPI process gets a unique seed

  for(unsigned int i=1; i<=N; ++i) {
    sample_rand(-512.,512.,dim,x);
    double f_i = secret_function(x[1],x[2]);
  //  printf( "Proccess % d, sample % i,  with f = % e \n", rank, i, f_i ) ;  
    if( f_i < min) {
      min = f_i;
    }
  }
  
// MPI_Barrier (MPI_COMM_WORLD);
//  printf("Process %d of %d local min = % e\n", rank, size, min );
 
 // double ring_time = omp_get_wtime();
  // ************ MPI MIN A Ring Topolog ************ 

//  double send_junk = min;
//  double rec_junk;
//  MPI_Status status;

//  if(rank==0) {
//    send_junk = min;
//    // only the process with rank ID = 0 will be in this block of code.
//    MPI_Send(&send_junk, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); //  send data to process 1
//    MPI_Recv(&rec_junk, 1, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD, &status); // receive data from process size-1
   

// //   printf("\n A Ring Topology Min = % e \n", rec_junk );
//    printf("\n%u, %i, %2.6e", N, size, rec_junk );
//  }
//  else if( rank == size-1) { 
//    MPI_Recv(&rec_junk, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status); // recieve data from process rank-1 (it "left" neighbor")
//    if(rec_junk < min){
//        send_junk = rec_junk;
//    }else{
//        send_junk = min;
//    }
//    MPI_Send(&send_junk, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send data to its "right neighbor", rank 0
//  }
//  else {
//    MPI_Recv(&rec_junk, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status); // recieve data from process rank-1 (it "left" neighbor")
//    if(rec_junk < min){
//        send_junk = rec_junk;
//    }else{
//        send_junk = min;
//    }
//    MPI_Send(&send_junk, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD); // send data to its "right neighbor" (rank+1)
//  }

  // ******************* MPI Min Reduction ******************* 
  double reduction_time = omp_get_wtime();

  MPI_Reduce(&min, &min_global, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  
  if( rank == 0) {
 //   printf("\n MPI Min Reduction Global Min = % e \n", min_global);
 //   printf(", %2.6e", min_global);
    printf("\n%u, %i, %2.6e", N, size,  min_global);
  }

  MPI_Finalize(); // programs should always perform a "graceful" shutdown

  double finish_time = omp_get_wtime();
 
  if( rank == 0 ) { 
//    printf("\n                   N : %i \n", N);
//    printf("                 size: %i \n", size );
//    printf(" Entire portion time : %f \n", (finish_time - start_time) );
//    printf(" Ring passing time   : %f \n", (reduction_time - ring_time) );
//    printf(" MPI Reduction time  : %f \n", (finish_time - reduction_time) );
//    printf(" Secret+Rnad functime: %f , (ring_time - start_time)/(size*N)  \n", (ring_time - start_time)/(size*N) ); 	 

//    printf(", %3.6f", (finish_time - start_time) );
//    printf(", %3.6f", (reduction_time - ring_time) );
//    printf(", %3.6f", (finish_time - reduction_time) );
//    printf(", %3.6f\n", (ring_time - start_time) );

    printf(", %3.6f", (finish_time - start_time) );
    printf(", %3.6f", (finish_time - reduction_time) );
    printf(", %3.6f\n", (reduction_time - start_time) );

//   printf(", %3.6f, %3.6f\n", (ring_time - start_time),  (ring_time - start_time)/(size*N) );

  }
  return 0;
}

