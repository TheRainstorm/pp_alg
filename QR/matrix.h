#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

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
matrix get_MM(matrix MATRIX1, matrix MATRIX2);
matrix get_zero_matrix(int m, int n);
matrix get_random_matrix(int m, int n);
matrix get_identity_matrix(int n);
matrix get_copy_matrix(matrix MATRIX1);

void copy_matrix(matrix Dst, matrix Src);
void transpose_matrix(matrix MATRIX1);

void free_matrix(matrix A);
void free_matrix_n(int n, ...);

static inline int get_matrix_size(matrix MATRIX1){
    return MATRIX1.m*MATRIX1.n*sizeof(double);
}

matrix read_matrix(const char *file){
    FILE *fp = fopen(file, "r");
    matrix MATRIX1;
    fscanf(fp, "%d %d", &MATRIX1.m, &MATRIX1.n);
    MATRIX1.data = (double *)malloc(sizeof(double) * MATRIX1.m * MATRIX1.n);
    if(MATRIX1.data==NULL){
        printf("Error: read matrix malloc failed\n");
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

//M3 = M1 x M2
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

//return  M1 x M2
matrix get_MM(matrix MATRIX1, matrix MATRIX2){
    matrix MATRIX3;
    MATRIX3.m = MATRIX1.m;
    MATRIX3.n = MATRIX2.n;
    MATRIX3.data = NULL;
    if(MATRIX1.n != MATRIX2.m){
        return MATRIX3;   //size不一致
    }
    MATRIX3.data = (double *)malloc(MATRIX3.m*MATRIX3.n*sizeof(double));

    MM(MATRIX1, MATRIX2, MATRIX3);
    return MATRIX3;
}

matrix get_random_matrix(int m, int n){
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

matrix get_zero_matrix(int m, int n){
    matrix MATRIX1;
    MATRIX1.m = m;
    MATRIX1.n = n;
    MATRIX1.data = (double *)malloc(m*n*sizeof(double));
    memset(MATRIX1.data, 0, m*n*sizeof(double));

    return MATRIX1;
}

matrix get_identity_matrix(int n){
    matrix MATRIX1;
    MATRIX1.m = MATRIX1.n = n;
    MATRIX1.data = (double *)malloc(n*n*sizeof(double));
    memset(MATRIX1.data, 0, n*n*sizeof(double));
    for(int i=0; i<n; i++)
        MATRIX1(i, i) = 1.0;

    return MATRIX1;
}

matrix get_copy_matrix(matrix MATRIX1){
    matrix MATRIX2 = MATRIX1;
    MATRIX2.data = (double *)malloc(get_matrix_size(MATRIX1));

    memcpy(MATRIX2.data, MATRIX1.data, get_matrix_size(MATRIX1));

    return MATRIX2;
}

void copy_matrix(matrix Dst, matrix Src){
    memcpy(Dst.data, Src.data, get_matrix_size(Src));
}

void transpose_matrix(matrix MATRIX1){
    double tmp;
    for(int i=0; i<MATRIX1.m; i++){
		for(int j=i+1; j<MATRIX1.n; j++){
            tmp = MATRIX1(i, j);
			MATRIX1(i, j) = MATRIX1(j, i);
            MATRIX1(j, i) = tmp;
		}
    }
}

void free_matrix(matrix A){
    free(A.data);
}
void free_matrix_n(int n, ...){
    va_list arg_ptr;
    va_start(arg_ptr, n);
    for(int i=0; i<n; i++){
        matrix tmp = va_arg(arg_ptr, matrix);
        free(tmp.data);
    }
    va_end(arg_ptr);
}