#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

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
matrix get_MM(matrix MATRIX1, matrix MATRIX2);
matrix init_matrix_random(int m, int n);


int main(){
    matrix MATRIX1;
    matrix MATRIX2;
    matrix MATRIX3;
    ////读取文件
    MATRIX1 = read_matrix("input/A.txt");
    MATRIX2 = read_matrix("input/B.txt");

    if(MATRIX1.n != MATRIX2.m){
        printf("Matrix A, B not fit!\n");
        exit(1);
    }

    int M = MATRIX1.m;
    int L = MATRIX1.n;
    int N = MATRIX2.n;

    printf("read A B matrix: %d %d %d\n", M, L, N);

    // //随机初始化A, B矩阵
    // M = L = N = 6;
    // MATRIX1 = init_matrix_random(M, L);
    // MATRIX2 = init_matrix_random(L, N);
    // print_matrix("A", MATRIX1);
    // print_matrix("B", MATRIX2);

    long t = -clock();
    MATRIX3 = get_MM(MATRIX1, MATRIX1);
    t += clock();
    printf("Matrix Mult finish!\n");
    printf("Test elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    free(MATRIX1.data);
    free(MATRIX2.data);
    free(MATRIX3.data);
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

matrix get_MM(matrix MATRIX1, matrix MATRIX2){
    matrix MATRIX3;
    MATRIX3.m = MATRIX1.m;
    MATRIX3.n = MATRIX2.n;
    MATRIX3.data = NULL;
    if(MATRIX1.n != MATRIX2.m){
        return MATRIX3;   //size不一致
    }
    MATRIX3.data = (double *)malloc(MATRIX3.m*MATRIX3.n*sizeof(double));

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
    return MATRIX3;
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