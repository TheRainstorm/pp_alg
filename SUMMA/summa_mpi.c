#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>

typedef struct matrix
{
    int m, n;
    double *data;
}matrix;

#define MATRIX1(i, j) MATRIX1.data[(i)*MATRIX1.n + j]
#define MATRIX2(i, j) MATRIX2.data[(i)*MATRIX2.n + j]
#define MATRIX3(i, j) MATRIX3.data[(i)*MATRIX3.n + j]

matrix read_matrix(const char *file);
void write_matrix(const char *file, matrix MATRIX1);
void print_matrix(const char *s, matrix MATRIX1);
int test_MM(matrix MATRIX1, matrix MATRIX2, matrix MATRIX3, double precise);
int MM(matrix MATRIX1, matrix MATRIX2, matrix MATRIX3);
matrix init_matrix_random(int m, int n);



int M, L, N;    //矩阵A, B的3个维度
int p, my_rank;
int pr, pc;     //处理器划分的行列
int ar, ac, br, bc;     //A, B划分的子块维度
int pi, pj;     //处理器所处的行, 列

#define A(i, j) A[(i)*L + j]
#define B(i, j) B[(i)*N + j]
#define C(i, j) C[(i)*N + j]
#define a(i, j) a[(i) * ac + j]
#define b(i, j) b[(i) * bc + j]
#define c(i, j) c[(i) * bc + j]
#define p(i, j) (i*pc + j)  //获得i行j列的处理器的编号

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //变量声明
    double *A, *B, *C;
    double *a, *b, *c, *a_col, *b_row;
    double start_time, end_time;
    int row_rank, col_rank;

    //初始化A，B矩阵
    if(my_rank==0){
        matrix MATRIX1;
        matrix MATRIX2;
        ////读取文件
        MATRIX1 = read_matrix("input/A.txt");
        MATRIX2 = read_matrix("input/B.txt");

        if(MATRIX1.n != MATRIX2.m){
            printf("Matrix A, B not fit!\n");
            exit(1);
        }

        M = MATRIX1.m;
        L = MATRIX1.n;
        N = MATRIX2.n;

        printf("read A B matrix: %d %d %d\n", M, L, N);

        // //随机初始化A, B矩阵
        // M = L = N = 6;
        // MATRIX1 = init_matrix_random(M, L);
        // MATRIX2 = init_matrix_random(L, N);
        // print_matrix("A", MATRIX1);
        // print_matrix("B", MATRIX2);


        A = MATRIX1.data;
        B = MATRIX2.data;

        C = (double *)malloc(sizeof(double)*M*N);
    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(my_rank==0){
        printf("M: %d L: %d N:%d\n", M, L, N);
    }

    start_time = MPI_Wtime();

    //分配处理器, 尽可能使pr, pc接近, 且p-pr*pc小
    // pr = sqrt(p);
    // pc = p/pr;

    // pr = p;
    // pc = 1;
    // while(pc < pr){
    //     pr >>= 1;
    //     pc = p / pr;
    // }

    int pr_save = p;
    pr = p;
    pc = 1;
    while(pc < pr){
        pc++;
        pr = p/pc;
        if(pr*pc==p){
            pr_save = pr;
        }
    }
    pr = pr_save; pc = p/pr_save;


    //划分矩阵
    ar = M / pr;
    ac = L / pc;
    br = L / pr;
    bc = N / pc;

    pi = my_rank / pc;
    pj = my_rank % pc;

    if(my_rank==0){
        printf("p %d pr %d pc %d\nar %d ac %d br %d bc %d\n", p, pr, pc, ar, ac, br, bc);

        if(M%pr || L%pc || L%pr || N%pc){
            printf("Error, Should statisify: pr | M, pc | N, pr,pc | L\n");
            exit(1);
        }
    }

    a = (double *)malloc(sizeof(double) * ar * ac);
    b = (double *)malloc(sizeof(double) * br * bc);
    c = (double *)malloc(sizeof(double) * ar * bc);
    a_col = (double *)malloc(sizeof(double) * ar);
    b_row = (double *)malloc(sizeof(double) * bc);

    if(a == NULL || b == NULL || c == NULL || a_col == NULL || b_row == NULL){
        printf("rank: %d malloc failed\n", my_rank);
    }

    //子通信域划分
    MPI_Comm row_comm, col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, my_rank/pc, my_rank%pc, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, my_rank%pc, my_rank/pc, &col_comm);
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);

    if(my_rank==0){
        printf("Distributing A, B!\n");
    }
    //A, B子块分发
    if(my_rank==0){
        for(int i=pr-1; i>=0; i--){ //倒序，因此，最后遍历到0号，进程，跳过MPI_Send，只保留复制块
            for(int j=pc-1; j>=0; j--){ 
                //复制块
                for(int k=0; k<ar; k++){
                    for(int l=0; l<ac; l++){
                        a(k, l) = A(i*ar + k, j*ac + l);
                    }
                }
                for(int k=0; k<ar; k++){
                    for(int l=0; l<ac; l++){
                        b(k, l) = B(i*br + k, j*bc + l);
                    }
                }
                if(!(i==0&&j==0)){
                    // printf("rank: %d send a: %d tag: %d\n", my_rank, p(i, j), 0);
                    MPI_Send(a, ar*ac, MPI_DOUBLE, p(i, j), 0, MPI_COMM_WORLD);
                    MPI_Send(b, br*bc, MPI_DOUBLE, p(i, j), 1, MPI_COMM_WORLD);
                }
            }
        }
    }else{
        MPI_Recv(a, ar*ac, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, br*bc, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    memset(c, 0, ar * bc*sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank==0){
        printf("SUMMA...!\n");
    }
    //主过程
    for(int k=0; k < L; k++){
        int col=k/ac, row=k/br; //，B的行

        if(row_rank == k/ac){   //位于k列所在A的列划分
            for(int i=0; i<ar; i++){
                a_col[i] = a(i, k%ac);
            }
            MPI_Bcast(a_col, ar, MPI_DOUBLE, row_rank, row_comm);
        }else{
            MPI_Bcast(a_col, ar, MPI_DOUBLE, k/ac, row_comm);
        }

        if(col_rank == k/br){   //位于k行所在B的行划分
            for(int j=0; j<bc; j++){
                b_row[j] = b(k%ar, j);
            }
            MPI_Bcast(b_row, bc, MPI_DOUBLE, col_rank, col_comm);
        }else{
            MPI_Bcast(b_row, bc, MPI_DOUBLE, k/br, col_comm);
        }

        //计算C子块
        for(int i=0; i<ar; i++){
            for(int j=0; j<bc; j++){
                c(i, j) += a_col[i]*b_row[j];
            }
        }
    }

	MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank==0){
        printf("Receving sub block C!\n");
    }
    //接收C子块
    if(my_rank==0){
        for(int i=0; i<pr; i++){
            for(int j=0; j<pc; j++){
                if(!(i==0 && j==0)){
                    MPI_Recv(c, ar*bc, MPI_DOUBLE, p(i, j), p(i, j), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // printf("rank: %d recv: %d tag: %d\n", my_rank, p(i, j), p(i, j));
                }
                for(int k=0; k<ar; k++){
                    for(int l=0; l<bc; l++){
                        C(i*ar + k, j*bc + l) = c(k, l);
                    }
                }
            }
        }
    }else{
        MPI_Send(c, ar*bc, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if(my_rank==0){
        matrix MATRIX1 = {M, M, A};
        matrix MATRIX2 = {M, M, B};
        matrix MATRIX3 = {M, M, C};
        
        printf("Matrix finish!\n");
        printf("elapsed time: %f\n", end_time-start_time);
        
        write_matrix("output/C.txt", MATRIX3);
        // print_matrix("C", MATRIX3);
        
        // if(test_MM(MATRIX1, MATRIX2, MATRIX3, 1e-6)){
        //     printf("Test success!\n");
        // }else{
        //     printf("Test failed!\n");
        // }

        free(A);
        free(B);
        free(C);
    }

    free(a);
    free(b);
    free(c);
    free(a_col);
    free(b_row);
    MPI_Finalize();

    return 0;
}

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
            fprintf(fp, "%6.3f ", MATRIX1(i, j));
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
            printf("%6.3f ", MATRIX1(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

//验证矩阵乘法M1xM2=M3
int test_MM(matrix MATRIX1, matrix MATRIX2, matrix MATRIX3, double precise){
    if(MATRIX1.n != MATRIX2.m || MATRIX1.m != MATRIX3.m || MATRIX2.n != MATRIX3.n){
        return 0;   //size不一致
    }
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX2.n; j++){
            double sum = 0;
            for(int k=0; k<MATRIX1.n; k++){
                sum += MATRIX1(i, k) * MATRIX2(k, j);
            }
            if(fabs(MATRIX3(i, j)-sum)>precise){
                return 0;
            }
        }
    }
    return 1;
}

int MM(matrix MATRIX1, matrix MATRIX2, matrix MATRIX3){
    if(MATRIX1.n != MATRIX2.m || MATRIX1.m != MATRIX3.m || MATRIX2.n != MATRIX3.n){
        return 0;   //size不一致
    }
    for(int i=0; i<MATRIX1.m; i++){
		for(int j=0; j<MATRIX2.n; j++){
			MATRIX3(i, j) = 0;
		}
        for(int k=0; k<MATRIX1.n; k++){
            for(int j=0; j<MATRIX2.n; j++){
                MATRIX3(i, j) += (MATRIX1(i, k) * MATRIX2(k, j));
            }
        }
    }
    return 1;
}

matrix init_matrix_random(int m, int n){
    matrix MATRIX1;
    MATRIX1.m = m;
    MATRIX1.n = n;
    MATRIX1.data = (double *)malloc(m*n*sizeof(double));
    for(int i=0; i<m; i++)
		for(int j=0; j<n; j++){
			MATRIX1(i, j) = rand()%100/100.0;
        }
    return MATRIX1;
}