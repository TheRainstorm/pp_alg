#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N    1000

double A[N][N], A_copy[N][N];
double Q[N][N], R[N][N];

void init_array_d2(double A[][N]);
void read_matrix(double A[][N], const char *file);
void write_matrix(const char *file, double A[][N]);
void print_array_d2(const char *s, double A[][N]);
void MM(double A[][N], double B[][N], double C[][N]);

void QR(double A[][N], double Q[][N], double R[][N]){
    //设置Q为单位矩阵
    memset(Q, 0, sizeof(double)*N*N);
    for(int i=0; i<N; i++){
        Q[i][i] = 1;
    }

    double sq, c, s;
    double aj[N], qj[N], ai[N], qi[N];
    for(int j=0; j<N; j++){
        for(int i=j+1; i<N; i++){
            sq = sqrt(A[j][j]*A[j][j] + A[i][j]*A[i][j]);
            c = A[j][j] / sq;
            s = A[i][j] / sq;

            //更新A和Q的i, j两行。把两个for循环合并后调bug调了1天。。。
            for(int k=0; k<N; k++){
                aj[k] = c*A[j][k] + s*A[i][k];
                qj[k] = c*Q[j][k] + s*Q[i][k];
                ai[k] = (-s)*A[j][k] + c*A[i][k];
                qi[k] = (-s)*Q[j][k] + c*Q[i][k];
            }
            for(int k=0; k<N; k++){
                A[j][k] = aj[k];
                Q[j][k] = qj[k];
                A[i][k] = ai[k];
                Q[i][k] = qi[k];
            }

            // for(int k=0; k<N; k++){
            //     aj[k] = c*A[j][k] + s*A[i][k];
            //     ai[k] = (-s)*A[j][k] + c*A[i][k];
            // }
            // for(int k=0; k<N; k++){
            //     A[j][k] = aj[k];
            //     A[i][k] = ai[k];
            // }
        }
    }

    memcpy(R, A, sizeof(double)*N*N);
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            double tmp = Q[i][j];
            Q[i][j] = Q[j][i];
            Q[j][i] = tmp;
        }
    }
}

//验证A*B = C
int valid_MM(const double A[][N], const double B[][N], double C[][N], double precise){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double sum = 0;
            for(int k=0; k<N; k++){
                sum += A[i][k] * B[k][j];
            }
            if(fabs(C[i][j]-sum)>precise){
                return 0;
            }
        }
    }
    return 1;
}

int main(){
    // srand(0);
    // init_array_d2(A);    //随机初始化
    read_matrix(A, "input/A.txt");

    // print_array_d2("A", A); 

    memcpy(A_copy, A, sizeof(double)*N*N);

    long t;
    t = -clock();
    QR(A, Q, R);
    t += clock();

    printf("Get Q, R!\n");
    printf("elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    write_matrix("output/Q_serial.txt", Q);
    write_matrix("output/R_serial.txt", R);
    // print_array_d2("Q", Q);
    // print_array_d2("R", R);

    // t = -clock();
    // if(valid_MM(Q, R, A_copy, 1e-3)){
    //     printf("Test success!\n");
    // }else{
    //     printf("Test failed! A can't QR decompose\n");
    // }
    // t += clock();
    // printf("Test elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    // MM(Q, R, A_copy);
    // print_array_d2("QR", A_copy);

    return 0;
}

void print_array_d2(const char *s, double A[][N]){
    printf("%s: \n", s);
    for(int i=0; i<N; i++){
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

void MM(double A[][N], double B[][N], double C[][N]){
    int i, j, k;
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            C[i][j] = 0;
        }
        for(k=0; k<N; k++){
            for(j=0; j<N; j++){
                C[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
}