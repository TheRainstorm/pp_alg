#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"

//N > 500
void calculate_serial(double A[][N], double C[][N], double D[][N]){
    for(int i=1; i<=100; i++){
        for(int j=1; j<=100; j++){
            if(j <= i + 6){
                A[3*i + 2*j][2*j] = C[i][j] * 2; //S1
                D[i][j] = A[i-j+6][i+j];         //S2
            }
        }
    }
}

void calculate_parallel(double A[][N], double C[][N], double D[][N]){
    #pragma omp parallel for
    for(int i=1; i<=100; i++){
        for(int j=1; j<=100; j++){
            if(j <= i + 6){
                A[3*i + 2*j][2*j] = C[i][j] * 2; //S1
                D[i][j] = A[i-j+6][i+j];         //S2
            }
        }
    }
}

double A[N][N], C[N][N], D[N][N];
double A2[N][N], C2[N][N], D2[N][N];

int main(){
    srand(0);

    init_array_d2_double(A);
    init_array_d2_double(C);
    init_array_d2_double(D);

    memcpy(A2, A, sizeof(A));
    memcpy(C2, C, sizeof(C));
    memcpy(D2, D, sizeof(D));

    long t;
    t = -clock();
    calculate_serial(A, C, D);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    t = -clock();
    calculate_parallel(A2, C2, D2);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    campare_d2_double("A", A, A2);
    campare_d2_double("D", D, D2);
    return 0;
}