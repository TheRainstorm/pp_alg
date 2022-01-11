#ifndef MY_MPI_API
#define MY_MPI_API
#include <mpi.h>

void my_Bcast_loop(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
void my_Bcast_fast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

void my_scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

void my_gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

void my_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm);

void my_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm);

#endif