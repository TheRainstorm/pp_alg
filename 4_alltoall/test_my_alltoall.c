#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "my_MPI_API.h"
#define N 1024
int isend[N], irecv[N];

int main(int argc, char *argv[]){
    int world_size, world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Test alltoall
    printf("before: rank = %d \t",world_rank);
    for(int i=0; i<world_size; i++)
        printf("%2d ", isend[i] = i + world_size * world_rank);
    printf( "\n");
    MPI_Barrier(MPI_COMM_WORLD);

    const int REPETE_TIMES=10000;
    double start_time, end_time;
    start_time = MPI_Wtime();
    for (int i = 0; i < REPETE_TIMES; i++){
        my_alltoall(isend, 1, MPI_INT, irecv, 1, MPI_INT, MPI_COMM_WORLD);
        // MPI_Alltoall(isend, 1, MPI_INT, irecv, 1, MPI_INT, MPI_COMM_WORLD);
    }
    end_time = MPI_Wtime();
    
    printf("after: rank = %d \t",world_rank);
    for(int i=0; i<world_size; i++)
        printf("%2d ", irecv[i]);
    printf( "\n");

    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0){
        printf("elapsed time: %.6f\n", (end_time-start_time));
    }

	MPI_Finalize();
    return 0;
}