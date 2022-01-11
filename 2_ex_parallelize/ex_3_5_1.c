#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"

void calculate_serial(double A[], double B[], double C[], double D[]){
    for(int i=1; i<N; i++){
        A[i] = A[i] + B[i-1];   //S1
        B[i] = C[i-1]*2;        //S2
        C[i] = 1/B[i];          //S3
        D[i] = C[i] * C[i];     //S4
    }
}

void calculate_parallel(double A[], double B[], double C[], double D[]){   
    for(int i=1; i<N; i++){
        B[i] = C[i-1]*2;        //S2
        C[i] = 1/B[i];          //S3
    }

	#pragma omp parallel for
    for(int i=1; i<N; i++){
        A[i] = A[i] + B[i-1];     //S1
    }
    
	#pragma omp parallel for
    for(int i=1; i<N; i++){
        D[i] = C[i] * C[i];     //S4
    }
}

double A[N], B[N], C[N], D[N];
double A2[N], B2[N], C2[N], D2[N];

int main(){
    srand(0);

    init_array_d1_double(A, N);
    init_array_d1_double(B, N);
    init_array_d1_double(C, N);
    init_array_d1_double(D, N);

    memcpy(A2, A, sizeof(A));
    memcpy(B2, B, sizeof(B));
    memcpy(C2, C, sizeof(C));
    memcpy(D2, D, sizeof(D));

    long t;
    t = -clock();
    calculate_serial(A, B, C, D);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    t = -clock();
    calculate_parallel(A2, B2, C2, D2);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    campare_d1_double("A", A, A2, N);
    campare_d1_double("B", B, B2, N);
    campare_d1_double("C", C, C2, N);
    campare_d1_double("D", D, D2, N);
    return 0;
}