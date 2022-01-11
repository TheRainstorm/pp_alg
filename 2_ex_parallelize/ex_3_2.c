#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"

#define Na 10

void calculate_serial(int A[][N], int B[], int C[][N], int X[], int Y[]){
    for(int i=0; i<N; i++){
        X[i] = Y[i] + 10; //S1
        for(int j=1; j<N-1; j++){   //j start from 1
            B[j] = A[j][Na]; //S2
            for(int k=0; k<N-1; k++){
                A[j+1][k] = B[j] + C[j][k]; //S3
            }
            Y[i+j] = A[j+1][Na]; //S4
        }
    }
}

void calculate_parallel(int A[][N], int B[], int C[][N], int X[], int Y[]){
    for(int j=1; j<N-1; j++){
        B[j] = A[j][Na];    //S2
        #pragma omp parallel for
        for(int k=0; k<N-1; k++){   //vector
            A[j+1][k] = B[j] + C[j][k];     //S3
        }
    }

    for(int i=0; i<N; i++)
        #pragma omp parallel for
        for(int j=1; j<N-1; j++)    //vector
            Y[i+j] = A[j+1][Na];    //S4

    #pragma omp parallel for
    for(int i=0; i<N; i++){
        X[i] = Y[i] + 10; //S1
    }
}

int A[N][N], B[N], C[N][N], X[N], Y[2*N];
int A2[N][N], B2[N], X2[N], Y2[2*N];

int main(){
    srand(0);

    init_array_d2(A);
    init_array_d1(B, N);
    init_array_d2(C);
    init_array_d1(X, N);
    init_array_d1(Y, 2*N);

    memcpy(A2, A, sizeof(A));
    memcpy(B2, B, sizeof(B));
    memcpy(X2, X, sizeof(X));
    memcpy(Y2, Y, sizeof(Y));

    long t;
    t = -clock();
    calculate_serial(A, B, C, X, Y);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    t = -clock();
    calculate_parallel(A2, B2, C, X2, Y2);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    campare_d2("A", A, A2);
    campare_d1("B", B, B2, N);
    campare_d1("X", X, X2, N);
    campare_d1("Y", Y, Y2, 2*N);
    return 0;
}