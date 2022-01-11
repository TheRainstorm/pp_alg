#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

int M;          //矩阵为M*M
int my_rank, p;

#define A(i, j) A[(i)*M + j]
#define a(i, j) a[(i)*M + j]

typedef struct matrix
{
    int m, n;
    double *data;
}matrix;

#define MATRIX1(i, j) MATRIX1.data[(i)*MATRIX1.n + j]
#define MATRIX2(i, j) MATRIX2.data[(i)*MATRIX2.n + j]
#define MATRIX3(i, j) MATRIX3.data[(i)*MATRIX3.n + j]

matrix read_matrix(const char *file){
    FILE *fp = fopen(file, "r");
    matrix MATRIX1;
    fscanf(fp, "%d %d", &MATRIX1.m, &MATRIX1.n);
    MATRIX1.data = (double *)malloc(sizeof(double) * MATRIX1.m * MATRIX1.n);
    if(MATRIX1.data==NULL){
        printf("Error: read matrix too big\n");
    }
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            fscanf(fp, "%lf", &MATRIX1(i, j));
        }
    }
    fclose(fp);
    return MATRIX1;
}

void write_matrix(const char *file, matrix MATRIX1){
    FILE *fp = fopen(file, "w+");
    fprintf(fp, "%d %d\n", MATRIX1.m, MATRIX1.n);
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            fprintf(fp, "%6.3f", MATRIX1(i, j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void print_matrix(const char *s, matrix MATRIX1){
    printf("%s: \n", s);
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            printf("%6.3f", MATRIX1(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

//验证矩阵乘法M1xM2=M3
int test_MM(matrix MATRIX1, matrix MATRIX2, matrix MATRIX3){
    if(MATRIX1.n != MATRIX2.m || MATRIX1.m != MATRIX3.m || MATRIX2.n != MATRIX3.n){
        return 0;   //size不一致
    }
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX2.n; j++){
            double sum = 0;
            for(int k=0; k<MATRIX1.n; k++){
                sum += MATRIX1(i, k) * MATRIX2(k, j);
            }
            if(fabs(MATRIX3(i, j)-sum)>1e-6){
                return 0;
            }
        }
    }
    return 1;
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double *A, *L, *U, *A_copy;
    double start_time, end_time;

    if(my_rank==0){
        //读取矩阵A
        matrix MATRIX1 = read_matrix("input/A.txt");
        if(MATRIX1.m != MATRIX1.n){
            printf("Input matrix should be square!\n");
            exit(0);
        }
        M = MATRIX1.m;
        A = MATRIX1.data;
        A_copy = (double *)malloc(sizeof(double)*M*M);
        memcpy(A_copy, A, sizeof(double)*M*M);
        // print_matrix("A", MATRIX1);
        printf("read matrix: %dx%d\n", M, M);

        //分配L, U空间
        L = (double *)malloc(sizeof(double)*M*M);
        U = (double *)malloc(sizeof(double)*M*M);
    }

    start_time = MPI_Wtime();
    //0号进程将M发送给其它进程
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int m = M / p;  //每个进程获得m行
    int r = M % p;  
    if(my_rank < r){ //非整除情况下，r个进程为m+1行，(p-r)个进程为m行（按交叉划分）
        m++;
    }
    // printf("rank: %d M: %d m: %d p: %d r: %d\n",my_rank, M, m, p, r);

    //每个进程分配m*M空间用于存储子矩阵
    double *a = (double *)malloc(sizeof(double)*m*M);
    
    //每个进程的发送和接收主行缓冲区
    double *buf = (double *)malloc(sizeof(double)*M);
    
    //0号进程将A的每一行，分配给p个进程（按交叉划分）
    if(my_rank==0){
        for(int i=0; i<M; i++){
            if(i%p ==0){    //排除掉0号进程自己的
                memcpy(&a(i/p, 0), &A(i, 0), M*sizeof(double));
            }else{
				MPI_Send(&A(i, 0), M, MPI_DOUBLE, i%p, i/p, MPI_COMM_WORLD);
            }
        }
    }else{
        for(int i=0; i<m; i++){
            MPI_Recv(&a(i, 0), M, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if(my_rank==0){
        printf("A distribute finish\n");
    }

    // if(my_rank==2){
    //     printf("rank: %d \n", my_rank);
    //     matrix tmp = {m, M, a};
    //     print_matrix("a", tmp);
    // }

    //主要逻辑：遍历M个主行，主行所在进程将主行发送给其它进程进行变换
    for(int k=0; k<M; k++){
        int block, proc;
        block = k/p;
        proc = k%p;

        if(my_rank == proc){
            memcpy(buf, &a(block, 0), M*sizeof(double));
        }
        //接收主行
        MPI_Bcast(buf, M, MPI_DOUBLE, proc, MPI_COMM_WORLD);

        int start_block = block;
        if(my_rank <= proc){     //小于my_rank，更新block+1后的行。大于则更新block后的行
            start_block++;
        }

        //更新进程所拥有的行
        for(int i=start_block; i<m; i++){
            a(i, k) = a(i, k) /buf[k];
            for(int j=k+1; j<M; j++){
                a(i, j) = a(i, j) - a(i, k)*buf[j];
            }
        }
    }

    if(my_rank==0){
        printf("A update finish\n");
    }

    //0号进程接收子行，构成变换后的A整体
    if(my_rank==0){
        for(int k=0; k<M; k++){
            int block, proc;
            block = k/p;
            proc = k%p;
            if(proc==my_rank){
                memcpy(&A(k, 0), &a(block, 0), M*sizeof(double));
            }else{
                MPI_Recv(&A(k, 0), M, MPI_DOUBLE, proc, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }else{
        for(int i=0; i<m; i++){
            MPI_Send(&a(i, 0), M, MPI_DOUBLE, 0, i*p+my_rank, MPI_COMM_WORLD);
        }
    }

    //0号进程通过A，计算L, U
    if(my_rank==0){
        matrix MATRIX1 = {M, M, L};
        matrix MATRIX2 = {M, M, U};
        matrix MATRIX3 = {M, M, A_copy};
        for(int i=0; i<M; i++){
            for(int j=0; j<M; j++){
                if(j<i){
                    MATRIX1(i, j) = A(i, j);
                    MATRIX2(i, j) = 0;
                }else{
                    MATRIX2(i, j) = A(i, j);
                    MATRIX1(i, j) = 0;
                }
            }
            MATRIX1(i, i) = 1;
        }
        end_time = MPI_Wtime();
        printf("Get L, U!\n");
        printf("elapsed time: %f\n", end_time-start_time);

        // print_matrix("L", MATRIX1);
        // print_matrix("U", MATRIX2);
        
        // write_matrix("output/L.txt", MATRIX1);
        // write_matrix("output/U.txt", MATRIX2);
        
        // start_time = MPI_Wtime();
        // if(test_MM(MATRIX1, MATRIX2, MATRIX3)){
        //     printf("Test success!\n");
        // }else{
        //     printf("Test failed! A is Singular Matrix\n");
        // }
        // end_time = MPI_Wtime();
        // printf("Test LU elapsed time: %f\n", end_time-start_time);

        free(A);
        free(A_copy);
        free(L);
        free(U);
    }

    free(a);
    free(buf);
    MPI_Finalize();

    return 0;
} 