#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 1000

double A[N][N], A_copy[N][N];
double L[N][N], U[N][N];

void LU(double A[][N], double L[][N], double U[][N]){
    // double *A_buffer = (double *)malloc(sizeof(A));
    for(int k=0; k<N; k++){
        for(int i=k+1; i<N; i++){
            A[i][k] = A[i][k]/A[k][k];      //保存当前行和第k行的比例
            for(int j=k+1; j<N; j++){
                A[i][j] = A[i][j] - A[i][k]*A[k][j];
            }
        }
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(j<i){
                L[i][j] = A[i][j];
                U[i][j] = 0;
            }else{
                U[i][j] = A[i][j];
                L[i][j] = 0;
            }
        }
        L[i][i] = 1;
    }
}

//验证A*B = C
int valid_MM(const double A[][N], const double B[][N], double C[][N]){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double sum = 0;
            for(int k=0; k<N; k++){
                sum += A[i][k] * B[k][j];
            }
            if(fabs(C[i][j]-sum)>1e-6){
                return 0;
            }
        }
    }
    return 1;
}

void print_array_d2(const char *s, double A[][N]){
    printf("%s: \n", s);
    for(int i=0; i<1; i++){
        for(int j=0; j<N; j++){
            printf("%6.3f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void init_array_d2(double A[][N]){
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            A[i][j] = rand()%100/100.0;
}

void read_matrix(double A[][N], const char *file){
    FILE *fp = fopen(file, "r");
    int m, n;
    fscanf(fp, "%d %d", &m, &n);
    if (m != n || m!=N){
        printf("The input matrix should be %d square!\n", N);
        exit(0);
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            fscanf(fp, "%lf", &A[i][j]);    //ele_t
    fclose(fp);
}

void write_matrix(const char *file, double A[][N]){
    FILE *fp = fopen(file, "w+");
    fprintf(fp, "%d %d\n", N, N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fprintf(fp, "%6.3f", A[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int main(){
    // srand(0);
    // init_array_d2(A);    //随机初始化

    read_matrix(A, "input/A.txt");
    // print_array_d2("A", A); 
    memcpy(A_copy, A, sizeof(double)*N*N);

    long t;
    t = -clock();
    LU(A, L, U);
    t += clock();

    printf("Get L, U!\n");
    printf("elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    write_matrix("output/L_serial.txt", L);
    write_matrix("output/U_serial.txt", U);
    // print_array_d2("L", L);
    // print_array_d2("U", U);

    t = -clock();
    if(valid_MM(L, U, A_copy)){
        printf("Test success!\n");
    }else{
        printf("Test failed! A is Singular Matrix\n");
    }
    t += clock();
    printf("Test LU elapsed time: %f\n", (double)t/CLOCKS_PER_SEC);

    return 0;
}