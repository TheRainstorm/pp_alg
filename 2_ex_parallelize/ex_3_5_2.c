#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"

void calculate_serial(double A[], double B[], double C[], double D[]){
    for(int i=0; i<N; i++){
        A[i] = B[i] + C[i];               //S1
        D[i] = (A[i] + A[N - 1 - i])/2;     //S2
    }
}

void calculate_parallel(double A[], double B[], double C[], double D[]){
    int i;
	#pragma omp parallel for private(i)
	for (i = 0; i <= (N - 1)/ 2; i++) {
		A[i] = B[i] + C[i];
		D[i] = (A[i] + A[N - 1 - i]) / 2;
	}
	#pragma omp parallel for private(i)
	for (i = (N - 1)/ 2 + 1; i < N; i++) {
		A[i] = B[i] + C[i];
		D[i] = (A[i] + A[N - 1 - i]) / 2;
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