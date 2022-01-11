#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"

#define Na 10

void calculate_serial(int A[][N]){  //N > 302
    for(int i=1; i<=100; i++){
        for(int j=1; j<=50; j++){
            A[3*i+2][2*j-1] = A[5*j][i+3] + 2;
        }
    }
}

void calculate_parallel(int A[][N]){
    for(int i=1; i<=100; i++){
        #pragma omp parallel for
        for(int j=1; j<=50; j++){
            A[3*i+2][2*j-1] = A[5*j][i+3] + 2;
        }
    }
}

int A[N][N];
int A2[N][N];

int main(){
    srand(0);

    init_array_d2(A);

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

    campare_d2("A", A, A2);
    return 0;
}