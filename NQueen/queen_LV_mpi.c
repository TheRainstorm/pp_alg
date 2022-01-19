#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define BT_STEP 20  //剩余使用回溯法的步数

int p, my_rank;
int n;  //棋盘大小
int *board;
long long run_time;    //记录每个进程运行次数

/*枚举信息标签*/
enum msg_tag
{
    SUCCESS_TAG,
    SOLUTIONS_TAG,
    TERMINATE_TAG,
    // REQUEST_TAG,
    // RESPONSE_TAG
};

void get_parameter(int argc, char *argv[], int *n);     //从命令行参数中获得n
static inline int rand_int(int min, int max);
int queen_lv(int *board, int n, int step);
void write_board(const char *file, int board[], int n);
void print_board(int board[], int n);

void Queen_LV_master(){
    /**
     * @brief master主要任务
     * 1. 接收来自其它进程的获取任务请求
     * 2.  
     * 
     */
    MPI_Status status;
    int success;
    int msg = 0, tag;
    int node;
    MPI_Request recv_handler;
    int flag;

    //接收其它进程的完成消息
    MPI_Irecv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, SUCCESS_TAG, MPI_COMM_WORLD, &recv_handler);
    
    while(1){
        // run_time ++;
        // int step = rand_int(n-BT_STEP, n);
        // success = queen_lv(board, n, step);
        // if(success){
        //     break;
        // }

        MPI_Test(&recv_handler, &flag, &status);
        if(flag){
            node = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            printf("recv result: %d\n", node);
            MPI_Recv(board, n, MPI_INT, node, SOLUTIONS_TAG, MPI_COMM_WORLD, &status);
            break;
        }
    }

    //给其它进程发送终止信号
    for(int i=1; i<p; i++){
        MPI_Send(&msg, 1, MPI_INT, i, TERMINATE_TAG, MPI_COMM_WORLD);
    }
    printf("rank %d end\n", my_rank);
}

void Queen_LV_slave(){
    MPI_Status status;
    int success;
    int msg = 0;
    MPI_Request recv_handler, send_handler;
    int flag;

    //接收主进程的终止信息
    MPI_Irecv(&msg, 1, MPI_INT, 0, TERMINATE_TAG, MPI_COMM_WORLD, &recv_handler);
    
    while(1){
        run_time ++;
        int step = n;
        // int step = rand_int(n-BT_STEP, n);
        success = queen_lv(board, n, step);
        if(success){
            //给master进程发送SUCCESS信息
            printf("send result: %d\n", my_rank);
            MPI_Isend(&msg, 1, MPI_INT, 0, SUCCESS_TAG, MPI_COMM_WORLD, &send_handler);
            //给master进程发送解
            MPI_Isend(board, n, MPI_INT, 0, SOLUTIONS_TAG, MPI_COMM_WORLD, &send_handler);
            break;
        }

        MPI_Test(&recv_handler, &flag, &status);
        if(flag){
            break;
        }
    }
    printf("rank %d end\n", my_rank);
}

int main(int argc, char *argv[]){
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand(my_rank*time(NULL));
    
    double start_time, end_time;
    long long run_time_total;
    if(my_rank==0){
        get_parameter(argc, argv, &n);
        printf("n=%d\n", n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);   //将棋盘大小广播给每个进程
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);   //将进程数广播给每个进程
    board = (int *)malloc(n*sizeof(n)); //给棋盘分配空间

    start_time = MPI_Wtime();
    if(my_rank==0)
        Queen_LV_master();
    else
        Queen_LV_slave();
    end_time = MPI_Wtime();

    MPI_Reduce(&run_time, &run_time_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if(my_rank==0){
        printf("elpased time: %.6f\n", end_time-start_time);
        printf("run times: %lld\n", run_time_total);
        write_board("board.txt", board, n);
    }
    
	MPI_Finalize();
    return 0;
}

void get_parameter(int argc, char *argv[], int *n_ptr){
    int n = 4;
    switch (argc)
    {
    case 2:
        n = atoll(argv[1]);
        break;
    case 1:
        break;
    default:
        printf("Usage: queen n");
        exit(1);
        break;
    }
    *n_ptr = n;
}

//判断第row行皇后放置在第col列是否与前面row-1行冲突
static inline int is_safe(int board[], int n, int row, int col){
    for(int i = 0; i < row; i++){
        if(board[i] == col || abs(board[i] - col) == row - i){
            return 0;
        }
    }
    return 1;
}

//返回[min, max]之间的随机数
static inline int rand_int(int min, int max){
    return min + rand() % (max - min + 1);
}

/**
 * 回溯法放置第row行的皇后
 * 将遍历的节点数目存储在node_traversed中。返回是否成功
 */
int backtrace(int board[], int n, int row, long long *node_traversed){
    if(row == n){
        return 1;
    }
    for(int col = 0; col < n; col++){
        if(is_safe(board, n, row, col)){
            board[row] = col;
            (*node_traversed)++;      //放好1个皇后，算是搜索了1个节点
            if(backtrace(board, n, row+1, node_traversed)){
                return 1;
            }
        }
    }
    return 0;
}

//有stepVegas参数的Las Vegas算法。前stepVegas行随机放置，之后行采用回溯法放置
int queen_lv(int *board, int n, int step){
    long long node_traversed = 0;
    for(int i=0; i<step; i++){ //放置第i行皇后
        int valid_place=0;
        int place;
        for(int j=0; j<n; j++){
            if(is_safe(board, n, i, j)){            //判断如果放置在j列的话，是否和前i-1行冲突
                valid_place++;
                if(rand_int(1, valid_place)==1){    //最终放置在每个可选位置的概率相同，均为1/valid_place
                    place = j;
                }
            }
        }
        if(valid_place==0){     //搜索失败，返回0
            return 0;
        }
        board[i] = place;
        (node_traversed)++;      //放好1个皇后，算是搜索了1个节点
    }
    return backtrace(board, n, step, &node_traversed);;
}

void print_board(int board[], int n){
    int i, j;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(board[i] == j){
                printf("Q ");
            }
            else{
                printf("- ");
            }
        }
        printf("\n");
    }
}

void write_board(const char *file, int board[], int n){
    FILE *fp = fopen(file, "w+");
    fprintf(fp, "board n: %d\n", n);
    
    int i, j;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(board[i] == j){
                fprintf(fp, "Q ");
            }
            else{
                fprintf(fp, "- ");
            }
        }
        fprintf(fp, "\n");
    }
}