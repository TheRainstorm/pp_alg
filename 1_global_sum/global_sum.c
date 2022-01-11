#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


//二叉树形求和
//type = int
//root = 0
//op = SUM
void my_Reduce_sum1(const int *sendbuf, int *recvbuf, MPI_Comm comm)
{
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    int send = *sendbuf, recv;
    int from, to;

    int i;
    //reduce
    i = 1;
    while(i<comm_size){
        if(comm_rank%i==0){
            if(comm_rank%(i<<1)==0){
                // printf("rank %d recv from %d\n", comm_rank, comm_rank^i);
                MPI_Recv(&recv, 1, MPI_INT, comm_rank^i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                send += recv;
            }else{
                MPI_Send(&send, 1, MPI_INT, comm_rank^i, 0, MPI_COMM_WORLD);
                // printf("rank %d send to %d\n", comm_rank, comm_rank^i);
            }
        }
        i <<= 1;
    }

    recv = send;    //for root
    //bcast
    while(i!=1){
        i >>= 1;
        if(comm_rank%i==0){
            if(comm_rank%(i<<1)==0){
                MPI_Send(&send, 1, MPI_INT, comm_rank^i, 0, MPI_COMM_WORLD);
            }else{
                MPI_Recv(&recv, 1, MPI_INT, comm_rank^i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    *recvbuf = recv;
}

//蝶式求和(butterfly)
//type = int
//root = 0
//op = SUM
void my_Reduce_sum2(const int *sendbuf, int *recvbuf, MPI_Comm comm)
{
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    int send = *sendbuf, recv;
    int from, to;

    int i = 1;
    while(i<comm_size){
        from = to = comm_rank ^ i;
        // printf("rank %d sendrecv to %d\n", comm_rank, comm_rank^i);
        MPI_Sendrecv(&send, 1, MPI_INT, to, 0,
                     &recv, 1, MPI_INT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        send += recv;
        i <<= 1;
    }

    *recvbuf = send;
}

int main(int argc, char *argv[]){
    int world_size, world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int select;
    switch (argc)
    {
    case 2:
        select = atoi(argv[1]);
        break;
    case 1:
        select = 0;
        break;
    default:
        if(world_rank==0)
            printf("Usage: global_sum 0|1|2\n");
        exit(0);
        break;
    }

    int num, sum;
    num = world_rank;

    double start_time, end_time;

    int REPETE_TIMES = 10;
    start_time = MPI_Wtime();
    for(int i=0; i<REPETE_TIMES; i++){
        switch (select)
        {
        case 0:
            MPI_Reduce(&num, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            break;
        case 1:
            my_Reduce_sum1(&num, &sum, MPI_COMM_WORLD);
            break;
        case 2:
            my_Reduce_sum2(&num, &sum, MPI_COMM_WORLD);
        default:
            break;
        }
    }
    end_time = MPI_Wtime();

    if(world_rank==0){
        printf("sum=%d\n", sum);
    }

    if(world_rank==0){
        printf("elapsed time: %.6f\n", (end_time-start_time)/REPETE_TIMES);
    }
    
	MPI_Finalize();
    return 0;
}