#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BT_STEP 20  //剩余使用回溯法的步数

//判断第row行皇后放置在第col列是否与前面row-1行冲突
static inline int is_safe(int board[], int n, int row, int col){
    for(int i = 0; i < row; i++){
        if(board[i] == col || abs(board[i] - col) == row - i){
            return 0;
        }
    }
    return 1;
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

//返回[min, max]之间的随机数
static inline int rand_int(int min, int max){
    return min + rand() % (max - min + 1);
}

/**
 * 有stepVegas参数的Las Vegas算法。前stepVegas行随机放置，之后行采用回溯法放置
 * 
 * 将遍历的节点数目存储在node_traversed中。返回是否成功
 */
int NQueen_LV_stepVegas(int board[], int n, int stepVegas, long long *node_traversed){
    *node_traversed = 0;
    if(stepVegas < 0 || stepVegas > n){
        printf("stepVegas参数错误！\n");
        exit(1);
    }
    for(int i=0; i<stepVegas; i++){ //放置第i行皇后
        int valid_place=0;
        int place;
        for(int j=0; j<n; j++){
            //判断如果放置在j列的话，是否和前i-1行冲突
            if(is_safe(board, n, i, j)){
                valid_place++;
                if(rand_int(1, valid_place)==1){    //最终放置在每个可选位置的概率相同，均为1/valid_place
                    place = j;
                }
            }
        }
        if(valid_place==0){
            return 0;
        }
        board[i] = place;
        (*node_traversed)++;      //放好1个皇后，算是搜索了1个节点
    }
    
    return backtrace(board, n, stepVegas, node_traversed);
}

void write_board(const char *file, int board[], int n);
void print_board(int board[], int n);

int main(int argc, char const *argv[]){
    srand(time(NULL));

    long long n = 4;
    switch (argc)
    {
    case 2:
        n = atoll(argv[1]);
        break;
    case 1:
        break;
    default:
        exit(1);
        break;
    }
    printf("n=%lld\n", n);

    int *board;
    board = (int *)malloc(n*sizeof(int));

    long long t = -clock();
    int is_success = 0;
    long long cnt = 0;
    long long node_traversed, node_traversed_total= 0;
    long long LV_step;
    while(!is_success){
        cnt ++;
        LV_step = n;
        // LV_step = n-BT_STEP;
        // LV_step = rand_int(n-BT_STEP, n);
        is_success = NQueen_LV_stepVegas(board, n, LV_step, &node_traversed);
        node_traversed_total += node_traversed;
    }
    t += clock();
    printf("elpased time: %.6f\n", (double)t/CLOCKS_PER_SEC);
    printf("LV repete time: %lld\n", cnt);
    printf("node traversed total: %lld\n", node_traversed_total);

    write_board("board.txt", board, n);
    // print_board(board, n);
    
    free(board);
    return 0;
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