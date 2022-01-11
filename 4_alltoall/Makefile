CC=mpicc
CFLAGS= -L. -lmy_mpi

all: test_my_alltoall test_my_allgather

test_my_alltoall: libmy_mpi.a test_my_alltoall.c
	$(CC) test_my_alltoall.c $(CFLAGS) -o $@

test_my_allgather: libmy_mpi.a test_my_allgather.c
	$(CC) test_my_allgather.c $(CFLAGS) -o $@

libmy_mpi.a: my_MPI_API.c
	$(CC) -c my_MPI_API.c -o $@

clean:
	rm libmy_mpi.a test_my_alltoall test_my_allgather
