#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util.h"

void calculate_serial(double A[]){
    for(int i=1; i<N; i++){
        A[i] =  A[i-1] + 1;
    }
}

void calculate_parallel(double A[]){
    #pragma omp parallel for
    for(int i=1; i<N; i++){
        A[i] =  A[i-1] + 1;
    }
}

int main(){
    srand(0);

    double A[N];
    double A2[N];

    init_array_d1_double(A, N);

    memcpy(A2, A, sizeof(A));

    long t;
    t = -clock();
    calculate_serial(A);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    t = -clock();
    calculate_parallel(A2);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);

    campare_d1_double("A", A, A2, N);
    return 0;
}