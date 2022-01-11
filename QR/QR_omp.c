#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "matrix.h"
#include <omp.h>

#define A(i,j) A[(i) * M + j]
#define Q(i,j) Q[(i) * M + j]
#define R(i,j) R[(i) * M + j]

void update(int i, int j, double *A, double *Q, double *R, double *aj, double *qj, int M){
    double sq, c, s;
    sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
    c = A(j, j) / sq;
    s = A(i, j) / sq;
    int k;
    for (k = 0; k < M; k++){
        aj[k] = c * A(j, k) + s * A(i, k);
        qj[k] = c * Q(j, k) + s * Q(i, k);
        A(i, k) = -s * A(j, k) + c * A(i, k);
        Q(i, k) = -s * Q(j, k) + c * Q(i, k);
    }

    for (k = 0; k < M; k++){
        A(j, k) = aj[k];
        Q(j, k) = qj[k];
    }
}

void print_f(int tid, int *f, int M){
    printf("tid: %d \tf: ", tid);
    for(int i=0; i<M; i++){
        printf("%d ", f[i]);
    }
    printf("\n");
}

void QR_omp(double *A, double *Q, double *R, int M){
	// omp_set_num_threads(1);

    int p = omp_get_max_threads();
    printf("QR omp: threads %d\n", p);
    int m, r;
    double start, end;

    if(p>M){
        p = M;
        omp_set_num_threads(p);
    }
    m = M / p;  //每个线程处理的行数	
    r = M % p;  //最后一个线程 m + r
    printf("QR omp: p: %d M: %d m: %d r: %d\n", p, M, m, r);

    // double *aj = (double *)malloc(sizeof(double) * M);
    // double *qj = (double *)malloc(sizeof(double) * M);
	int *f = (int *)calloc(M, sizeof(int)); //column finish flag, set to zero
	// int *f = (int *)malloc(sizeof(int) * M); //column finish flag, set to zero
    // memset(f, 0, sizeof(int)*M);
    // int f[M];

	start = omp_get_wtime();
    #pragma omp parallel
    {
        int i, j, k, bias;
        double aj[M];
        double qj[M];
        int tid = omp_get_thread_num();

        bias = tid==(p-1) ? r : 0;
        for(j=0; j<tid*m; j++){     //矩形
            while(f[j]!=tid){
                // printf("tid:%d wait\n", tid);
            }

            for(i=tid*m; i<(tid+1)*m + bias; i++){
                update(i, j, A, Q, R, aj, qj, M);
            }

            f[j]=tid + 1;
            // printf("tid:%d Rec: j:%d, f[j]:%d\n", tid, j, f[j]);
            // print_f(tid, f, M);
        }

        for(j = tid*m; j<(tid+1)*m-1+bias; j++){ //三角
            for(i=j+1; i<(tid+1)*m+bias; i++){
                update(i, j, A, Q, R, aj, qj, M);
            }
            f[j] = tid + 1;
            // printf("tid:%d Tri: j:%d, f[j]:%d\n", tid, j, f[j]);
            // print_f(tid, f, M);
        }
        f[(tid+1)*m-1] = tid + 1;
        // print_f(tid, f, M);


		#pragma omp barrier
        // printf("tid:%d finish\n", tid);
        // if(tid==0){
        //     matrix mA={M, M, A};
        //     matrix mQ={M, M, Q};
        //     matrix mR={M, M, R};
        //     transpose_matrix(mQ);
        //     copy_matrix(mR, mA);
        //     if(test_MM(mQ, mR, A_copy, 1e-3)){
        //         printf("QR Test Success!\n");
        //     }else{
        //         printf("QR Test Faild!\n");
        //     }
        // }
    }
    // printf("hello\n");

    matrix mA={M, M, A};
    matrix mQ={M, M, Q};
    matrix mR={M, M, R};
    transpose_matrix(mQ);
    copy_matrix(mR, mA);
    // print_matrix("A", mA);
    // print_matrix("Q", mQ);
    // print_matrix("R", mR);
    
    // free(aj);
    // free(qj);
    // free(f);
}

int main(int argc, char **argv)
{
    //设置A数组
    matrix A = read_matrix("input/A.txt");
    if(A.m != A.n){
        printf("Input matrix should be square!\n");
        exit(0);
    }
    // print_matrix("A", A);
    printf("read matrix: %dx%d\n", A.m, A.m);

    matrix Q = get_identity_matrix(A.m);
    matrix R = get_zero_matrix(A.m, A.m);
    matrix A_copy = get_copy_matrix(A);
    
    // long t;
    // t = -clock();
    double start, end;
	start = omp_get_wtime();
    QR_omp(A.data, Q.data, R.data, A.m);
	end = omp_get_wtime();
    printf("elapsed time: %f\n", end-start);

    printf("get Q,R!\n");
    // if(test_MM(Q, R, A_copy, 1e-3)){
    //     printf("QR Test Success!\n");
    // }else{
    //     printf("QR Test Faild!\n");
    // }

    write_matrix("output/Q.txt", Q);
    write_matrix("output/R.txt", R);

    free_matrix_n(4, A, Q, R, A_copy);
	return(0);
}