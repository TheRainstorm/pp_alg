#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

unsigned long hash(char *str)
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

int main(int argc, char *argv[]){
    int world_size, world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    long long msg_size;
    int bcast_select;
    msg_size = 10000;
    bcast_select = 2;
    switch (argc)
    {
    case 3:
        bcast_select = atoi(argv[2]);
    case 2:
        msg_size = atoll(argv[1]);
        break;
    case 1:
        break;
    default:
        if(world_rank==0)
            printf("Usage: Bcast 0|1|2 [msg size]\n");
        exit(0);
        break;
    }

    /**
     * 使用节点主机名划分子通信域
     */
    MPI_Comm node_comm;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
    ////简单起见，暂时不使用node name
    // int color, key;
    // color = world_rank%4;
    // key = world_rank/4;
    // MPI_Comm_split(MPI_COMM_WORLD, color, key, &node_comm);
    MPI_Comm_split(MPI_COMM_WORLD, hash(processor_name), world_rank, &node_comm);   //?直接使用world_rank貌似没问题
    int node_size, node_rank;
    MPI_Comm_size(node_comm, &node_size);
	MPI_Comm_rank(node_comm, &node_rank);

    /**
     * 将节点通信域rank=0的进程，组成head通信域
     */
    int *node_ranks = (int *)malloc(sizeof(int)*world_size);
    int *head_ranks = (int *)malloc(sizeof(int)*world_size);
    MPI_Allgather(&node_rank, 1, MPI_INT, node_ranks, 1, MPI_INT, MPI_COMM_WORLD);  //每个进程都获得其他进程的节点通信域rank

    int idx=0;
    for(int i=0; i<world_size; i++){    //找到node_rank为0对应的world_rank
        if(node_ranks[i]==0){
            head_ranks[idx++] = i;
        }
    }

    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group head_group;
    MPI_Group_incl(world_group, idx, head_ranks, &head_group);  //构造head group
    MPI_Comm head_comm;
    MPI_Comm_create_group(MPI_COMM_WORLD, head_group, 0, &head_comm);

    int head_size=-1, head_rank=-1;
    if(MPI_COMM_NULL != head_comm){
        MPI_Comm_size(head_comm, &head_size);
        MPI_Comm_rank(head_comm, &head_rank);
    }
    
    /**
     * bcast
     */
    char *send_buf = (char *)malloc(msg_size);
    if(world_rank==0){
        sprintf(send_buf, "hello from root\n");
    }else{
        sprintf(send_buf, "NOT RECEIVED!!!\n");
    }

    double start_time, end_time;
    if(MPI_COMM_NULL != head_comm){ //head通信域中的广播
        switch (bcast_select){
            case 0:
                MPI_Bcast(send_buf, msg_size, MPI_CHAR, 0, head_comm);
                break;
            case 1:
                my_Bcast_loop(send_buf, msg_size, MPI_CHAR, 0, head_comm);
                break;
            case 2:
                my_Bcast_fast(send_buf, msg_size, MPI_CHAR, 0, head_comm);
                break;
            default:
                break;
        }
    }
    start_time = MPI_Wtime();
    int REPETE_TIMES = 1000;
    for(int i=0; i<REPETE_TIMES; i++){
        //在每个node通信域中的广播
        switch (bcast_select){
            case 0:
                MPI_Bcast(send_buf, msg_size, MPI_CHAR, 0, node_comm);
                break;
            case 1:
                my_Bcast_loop(send_buf, msg_size, MPI_CHAR, 0, node_comm);
                break;
            case 2:
                my_Bcast_fast(send_buf, msg_size, MPI_CHAR, 0, node_comm);
                break;
            default:
                break;
        }
    }
    end_time = MPI_Wtime();

    printf("Node Name: %s, WORLD RANK/SIZE: %2d/%2d --- NODE RANK/SIZE: %2d/%2d --- HEAD RANK/SIZE: %2d/%2d Rcv Msg: %s",
        processor_name, world_rank, world_size, node_rank, node_size, head_rank, head_size, send_buf);

    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0){
        printf("elapsed time: %.6f\n", (end_time-start_time)/REPETE_TIMES);
    }

    MPI_Finalize();
    return 0;
}