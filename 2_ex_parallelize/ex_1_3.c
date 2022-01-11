#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util.h"

void calculate_serial(double A[], double B[], double C[], double D[]){
    for(int i=0; i<N-1; i++){
        A[i] = B[i] + C[i+1];
        C[i] = A[i]*D[i];
    }
}

void calculate_parallel(double A[], double B[], double C[], double D[]){
    #pragma omp parallel for
    for(int i=0; i<N-1; i++){
        A[i] = B[i] + C[i+1];
        C[i] = A[i]*D[i];
    }
}

int main(){
    srand(0);

    double A[N], B[N], C[N], D[N];
    double A2[N], C2[N];

    init_array_d1_double(A, N);
    init_array_d1_double(B, N);
    init_array_d1_double(C, N);
    init_array_d1_double(D, N);

    memcpy(A2, A, sizeof(A));
    memcpy(C2, C, sizeof(C));

    long t;
    t = -clock();
    calculate_serial(A, B, C, D);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    t = -clock();
    calculate_parallel(A2, B, C2, D);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    campare_d1_double("A", A, A2, N);
    campare_d1_double("C", C, C2, N);
    return 0;
}