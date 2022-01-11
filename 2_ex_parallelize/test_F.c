#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 10

void init_array_d1(int A[], int n);
void init_array_d2(int A[][N]);
void campare_d1(const char *s, int A[], int B[], int n);
void campare_d2(const char *s, int A[][N], int B[][N]);

void init_array_d1_double(double A[], int n);
void init_array_d2_double(double A[][N]);
void campare_d1_double(const char *s, double A[], double B[], int n);
void campare_d2_double(const char *s, double A[][N], double B[][N]);

void print_array_d1(const char *s, int A[], int n);

void calculate_serial(int F[]){
    for(int i=0; i<N-4; i++){
        F[i] = F[i+1] + F[i+2] + F[i+3] + F[i+4];
    }
}

void calculate_parallel(int F[], int buf[]){
    #pragma omp parallel for
    for(int i=0; i<N-4; i++){
        buf[i] = F[i+1] + F[i+2] + F[i+3] + F[i+4];
    }
    #pragma omp parallel for
    for(int i=0; i<N-4; i++){
        F[i] = buf[i];
    }
}

int F[N];
int F2[N];
int buf[N];

int main(){
    srand(0);

    init_array_d1(F, N);
    print_array_d1("F", F, N);

    memcpy(F2, F, sizeof(F));

    long t;
    t = -clock();
    calculate_serial(F);
    t += clock();
    printf("Serial: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);
    print_array_d1("F", F, N);

    t = -clock();
    calculate_parallel(F2, buf);
    t += clock();
    printf("Parallel: elapsed time = %f\n", (double)t/CLOCKS_PER_SEC);
    print_array_d1("F2", F2, N);

    campare_d1("F", F, F2, N);
    return 0;
}


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
