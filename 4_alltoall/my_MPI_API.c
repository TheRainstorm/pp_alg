#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "my_MPI_API.h"

static unsigned long hash(char *str)
{
//https://stackoverflow.com/a/7666577/9933066
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

void my_Bcast_loop(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    if(root == comm_rank){
        for(int i=0; i<comm_size; i++){
            if(i!=root){
                MPI_Send(buffer, count, datatype, i, 0, comm);
            }
        }
    }else{
        MPI_Recv(buffer, count, datatype, root, 0, comm, MPI_STATUS_IGNORE);
    }
}

void my_Bcast_fast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    //使用ranks数组，对rank进行间接访问，从而对root非0的情况进行统一
    int *ranks = (int *)malloc(comm_size*sizeof(int));
    for(int i=0; i<comm_size; i++){
        ranks[i] = i;
    }
    if(root!=0){
        ranks[0] = root;
        ranks[root] = 0;
    }

    //某一轮需要在length个进程间发送消息，此时：
    //p —> p + length/2
    long long length = 2, offset;
    while(length<=comm_size){
        offset = length>>1;
        for(int i=offset; i<length && i<comm_size; i++){
            if(ranks[i-offset]==comm_rank){
                MPI_Send(buffer, count, datatype, ranks[i], 0, comm);
            }
            if(ranks[i]==comm_rank){
                MPI_Recv(buffer, count, datatype, ranks[i-offset], 0, comm, MPI_STATUS_IGNORE);
            }
        }
        length <<= 1;
    }
}

void my_scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    //root send: O(N)
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);
    int stype_size, rtype_size;
    MPI_Type_size(sendtype, &stype_size);
    MPI_Type_size(recvtype, &rtype_size);

    if(comm_rank==root){
        for(int i=0; i<comm_size; i++){     //send root -> i
            if(i==root)
                memcpy(recvbuf, sendbuf+i*sendcount*stype_size, stype_size*sendcount);
            else
                MPI_Send(sendbuf+i*sendcount*stype_size, sendcount, sendtype, i, 0, comm);
        }
    }else{      //recv root -> comm_rank
        MPI_Recv(recvbuf, recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
    }
}

void my_gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    //root recv: O(N)
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);
    int stype_size, rtype_size;
    MPI_Type_size(sendtype, &stype_size);
    MPI_Type_size(recvtype, &rtype_size);

    if(comm_rank==root){
        for(int i=0; i<comm_size; i++){     //send i -> root
            if(i==root)
                memcpy(recvbuf+i*recvcount*rtype_size, sendbuf, stype_size*sendcount);
            else
                MPI_Recv(recvbuf+i*recvcount*rtype_size, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }else{      //send comm_rank -> root
        MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
    }
}


// void my_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
//                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
// {
//     //不会死锁
//     int comm_size, comm_rank;
//     MPI_Comm_size(comm, &comm_size);
//     MPI_Comm_rank(comm, &comm_rank);

//     for(int i=0; i<comm_size; i++){         //gather i
//         my_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, i, comm);
//     }
// }

void my_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    //should be O(N)
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);
    int stype_size, rtype_size;
    MPI_Type_size(sendtype, &stype_size);
    MPI_Type_size(recvtype, &rtype_size);

    for(int i=0; i<comm_size; i++){         //gather i (i.e. root=i)
        if(comm_rank==i){   //this process is root
            for(int j=0; j<comm_size; j++){     //send j -> root
                if(i==j)
                    memcpy(recvbuf+j*recvcount*rtype_size, sendbuf, stype_size*sendcount);
                else
                    //store sender j's data at recvbuf[j]                        //from sender j
                    MPI_Recv(recvbuf+j*recvcount*rtype_size, recvcount, recvtype, j, 0, comm, MPI_STATUS_IGNORE);
            }
        }else{      //send comm_rank -> root
            MPI_Send(sendbuf, sendcount, sendtype, i, 0, comm);
        }
    }
}

// void my_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
//     int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
// {
//     int comm_size, comm_rank;
//     MPI_Comm_size(comm, &comm_size);
//     MPI_Comm_rank(comm, &comm_rank);
//     int rtype_size;
//     MPI_Type_size(sendtype, &rtype_size);

//     for(int i=0; i<comm_size; i++){         //i scatter
//         my_scatter(sendbuf, sendcount, sendtype, recvbuf+i*recvcount*rtype_size, recvcount, recvtype, i, comm);
//     }
// }

void my_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);
    int stype_size, rtype_size;
    MPI_Type_size(sendtype, &stype_size);
    MPI_Type_size(recvtype, &rtype_size);

    for(int i=0; i<comm_size; i++){         //i scatter (root=i)
        if(comm_rank==i){
            for(int j=0; j<comm_size; j++){     //send root -> i
                if(i==j)
                    memcpy(recvbuf+i*recvcount*rtype_size, sendbuf+j*sendcount*stype_size, stype_size*sendcount);
                else
                    MPI_Send(sendbuf+j*sendcount*stype_size, sendcount, sendtype, j, 0, comm);
            }
        }else{      //recv root -> comm_rank
            MPI_Recv(recvbuf+i*recvcount*rtype_size, recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }
}
