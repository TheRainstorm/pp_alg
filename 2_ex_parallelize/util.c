#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void init_array_d1(int A[], int n){
    for(int i=0; i<n; i++){
        A[i] = rand()%100;
    }
}

void init_array_d2(int A[][N]){
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            A[i][j] = rand()%100;
}

void campare_d1(const char *s, int A[], int B[], int n){
    for(int i=0; i<n; i++){
        if(A[i]!=B[i]){
            printf("%s: Diff: %d %d %d\n",s, i, A[i], B[i]);
            return;
        }
    }
    printf("%s: Correct\n", s);
}

void campare_d2(const char *s, int A[][N], int B[][N]){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(A[i][j]!=B[i][j]){
                printf("%s: Diff: %d %d %d %d\n", s, i, j, A[i][j], B[i][j]);
                return;
            }
        }
    }
    printf("%s: Correct\n", s);
}


void init_array_d1_double(double A[], int n){
    for(int i=0; i<n; i++){
        A[i] = rand()%100/100.0;
    }
}

void init_array_d2_double(double A[][N]){
    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            A[i][j] = rand()%100/100.0;
}

void campare_d1_double(const char *s, double A[], double B[], int n){
    for(int i=0; i<n; i++){
        if(fabs(A[i]- B[i])>1e-6){
            printf("%s: Diff: %d %f %f\n", s, i, A[i], B[i]);
            return;
        }
    }
    printf("%s: Correct\n", s);
}

void campare_d2_double(const char *s, double A[][N], double B[][N]){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(fabs(A[i][j]- B[i][j])>1e-6){
                printf("%s: Diff: %d %d %f %f\n", s, i, j, A[i][j], B[i][j]);
                return;
            }
        }
    }
    printf("%s: Correct\n", s);
}

void print_array_d1(const char *s, int A[], int n){
    printf("%s\n", s);
    for(int i=0; i<N; i++){
        printf("%d ", A[i]);
    }
    printf("\n");
}