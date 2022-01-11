#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#define A(x,y) A[(x) * M + y]
#define A0(x,y) A0[(x) * M + y]
#define Q(x,y) Q[(x) * M + y]
#define R(x,y) R[(x) * M + y]

#define THREADS_NUM 24
#define CLK_TCK 1000

double temp;
double *A, *A0;
double *ai, *qi, *aj, *qj;
double *R;
double *Q;
int p, m, r, M;

typedef struct matrix
{
    int m, n;
    double *data;
}matrix;

#define MATRIX1(i, j) MATRIX1.data[(i)*MATRIX1.n + j]
#define MATRIX2(i, j) MATRIX2.data[(i)*MATRIX2.n + j]
#define MATRIX3(i, j) MATRIX3.data[(i)*MATRIX3.n + j]

matrix read_matrix(const char *file){
    FILE *fp = fopen(file, "r");
    matrix MATRIX1;
    fscanf(fp, "%d %d", &MATRIX1.m, &MATRIX1.n);
    MATRIX1.data = (double *)malloc(sizeof(double) * MATRIX1.m * MATRIX1.n);
    if(MATRIX1.data==NULL){
        printf("Error: read matrix too big\n");
    }
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            fscanf(fp, "%lf", &MATRIX1(i, j));
        }
    }
    fclose(fp);
    return MATRIX1;
}

void write_matrix(const char *file, matrix MATRIX1){
    FILE *fp = fopen(file, "w+");
    fprintf(fp, "%d %d\n", MATRIX1.m, MATRIX1.n);
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            fprintf(fp, "%6.3f", MATRIX1(i, j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void print_matrix(const char *s, matrix MATRIX1){
    printf("%s: \n", s);
    for(int i=0; i<MATRIX1.m; i++){
        for(int j=0; j<MATRIX1.n; j++){
            printf("%6.3f", MATRIX1(i, j));
        }
        printf("\n");
    }
    printf("\n");
}
void RandomMatrix(double *A, double *A0, int n) {
	
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			A(i, j) = rand() % 1000;
			A0(i, j) = A(i, j);
		}
}

int TestQR(double *A0, double *Q, double *R, int n) {
	double tmp;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			tmp = 0;
			for (int k = 0; k < n; k++)
				tmp += Q(i, k) * R(k, j);
			if (fabs(A0(i, j) - tmp) > 0.01)
				return 0;
		}
	return 1;
}


int main(int argc, char **argv)
{
	int z;
	int i, j, k, nthreads, tid;
	double c, s, sq;
	int *f;
	double start, end;
	double ParTime;

	
	omp_set_num_threads(THREADS_NUM);
	// srand(0);

    matrix MATRIX1 = read_matrix("input/A.txt");
    if(MATRIX1.m != MATRIX1.n){
    printf("Input matrix should be square!\n");
    exit(0);
    }
    M = MATRIX1.m;
    A = MATRIX1.data;
    A0 = (double *)malloc(sizeof(double)*M*M);
    memcpy(A0, A, sizeof(double)*M*M);
    // print_matrix("A", MATRIX1);
    printf("read matrix: %dx%d\n", M, M);

    Q = (double *)malloc(sizeof(double) * M * M);
	R = (double *)malloc(sizeof(double) * M * M);

	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			if (i == j)
				Q(i, j) = 1.0;
			else
				Q(i, j) = 0.0;

	qi = (double *)malloc(sizeof(double) * M);
	qj = (double *)malloc(sizeof(double) * M);
	aj = (double *)malloc(sizeof(double) * M);
	ai = (double *)malloc(sizeof(double) * M);
	f = (int *)malloc(sizeof(int) * M);

	if (f == NULL || qi == NULL || qj == NULL || ai == NULL || aj == NULL)
		printf("memory allocation is wrong\n");

	for (i = 0; i < M; i++) {
		f[i] = 0;
	}

	start = omp_get_wtime();
	
	// 并行
	if (THREADS_NUM > 1) {
		#pragma omp parallel firstprivate(tid, i, j, k, s, sq, c, ai, aj, qi, qj) shared(A, Q, R, f, m, p, r, nthreads)
		{
			tid = omp_get_thread_num();
			if (tid == 0) {
				nthreads = omp_get_num_threads();
				p = nthreads;
				m = M / p;	// 每个线程处理的任务数	
				r = 0;	// 如果不能均分，最后一个线程处理的任务数
				if (M % p != 0) {
					m++;
					if ((M / m) < p) {	// 不需要用到这么多线程
						p = M / m;
						if (M % m != 0) {
							p++;
							r = M % m;
						}
					}
				}
			}
			#pragma omp barrier

			if (tid == 0) {
				for (j = 0; j < m - 1; j++)
				{
					for (i = j + 1; i < m; i++)
					{
						sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
						c = A(j, j) / sq;  s = A(i, j) / sq;

						for (k = 0; k < M; k++)
						{
							aj[k] = c * A(j, k) + s * A(i, k);
							qj[k] = c * Q(j, k) + s * Q(i, k);
							ai[k] = -s * A(j, k) + c * A(i, k);
							qi[k] = -s * Q(j, k) + c * Q(i, k);
						}

						for (k = 0; k < M; k++)
						{
							A(j, k) = aj[k];
							Q(j, k) = qj[k];
							A(i, k) = ai[k];
							Q(i, k) = qi[k];
						}
					}
					f[j] = 1;
					//printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
				}
				f[m - 1] = 1;	// 本处理器最后一行
				//printf("tid:%d, f:%d, j:%d\n", tid, f[m-1], m-1);
			}

			if (tid < p - 1 && tid > 0)
			{
				for (j = 0; j < tid * m; j++)	// 接收之前的主行
				{
					while (f[j] != tid);	// 等待前面的主行计算完成

					for (i = tid * m; i < (tid + 1) * m; i++)
					{
						sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
						c = A(j, j) / sq;  s = A(i, j) / sq;

						for (k = 0; k < M; k++)
						{
							aj[k] = c * A(j, k) + s * A(i, k);
							qj[k] = c * Q(j, k) + s * Q(i, k);
							ai[k] = -s * A(j, k) + c * A(i, k);
							qi[k] = -s * Q(j, k) + c * Q(i, k);
						}

						for (k = 0; k < M; k++)
						{
							A(j, k) = aj[k];
							Q(j, k) = qj[k];
							A(i, k) = ai[k];
							Q(i, k) = qi[k];
						}
					}
					//printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
					f[j] = tid + 1;
					// printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
				}

				for (j = tid * m; j < (tid + 1) * m - 1; j++)
				{

					for (i = j + 1; i < (tid + 1) * m; i++)
					{
						sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
						c = A(j, j) / sq;
						s = A(i, j) / sq;

						for (k = 0; k < M; k++)
						{
							aj[k] = c * A(j, k) + s * A(i, k);
							qj[k] = c * Q(j, k) + s * Q(i, k);
							ai[k] = -s * A(j, k) + c * A(i, k);
							qi[k] = -s * Q(j, k) + c * Q(i, k);
						}

						for (k = 0; k < M; k++)
						{
							A(j, k) = aj[k];
							Q(j, k) = qj[k];
							A(i, k) = ai[k];
							Q(i, k) = qi[k];
						}
					}
					//printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
					f[j] = tid + 1;
					// printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
				}
				//printf("tid:%d, f:%d, j:%d\n", tid, f[(tid + 1)*m - 1], (tid + 1)*m - 1);
				f[(tid + 1) * m - 1] = tid + 1;
				// printf("tid:%d, f:%d, j:%d\n", tid, f[(tid + 1)*m - 1], (tid+1)*m-1);
			}

			if (tid == p - 1) {	// p-1
				for (j = 0; j < tid * m; j++)
				{
					while (f[j] != p - 1);

					for (i = tid * m; i < (tid + 1) * m; i++)
					{
						if (!(r != 0 && i - tid * m >= r)) {
							sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
							c = A(j, j) / sq;  s = A(i, j) / sq;

							for (k = 0; k < M; k++)
							{
								aj[k] = c * A(j, k) + s * A(i, k);
								qj[k] = c * Q(j, k) + s * Q(i, k);
								ai[k] = -s * A(j, k) + c * A(i, k);
								qi[k] = -s * Q(j, k) + c * Q(i, k);
							}

							for (k = 0; k < M; k++)
							{
								A(j, k) = aj[k];
								Q(j, k) = qj[k];
								A(i, k) = ai[k];
								Q(i, k) = qi[k];
							}
						}
					}
					f[j] = tid + 1;
					// printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
				}

				for (j = tid * m; j < (tid + 1) * m - 1; j++)
				{

					if (!(r != 0 && j - tid * m >= r - 1)) {
						for (i = j + 1; i < (tid + 1) * m; i++)
						{
							if (!(r != 0 && i - tid * m >= r)) {
								sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
								c = A(j, j) / sq;
								s = A(i, j) / sq;

								for (k = 0; k < M; k++)
								{
									aj[k] = c * A(j, k) + s * A(i, k);
									qj[k] = c * Q(j, k) + s * Q(i, k);
									ai[k] = -s * A(j, k) + c * A(i, k);
									qi[k] = -s * Q(j, k) + c * Q(i, k);
								}

								for (k = 0; k < M; k++)
								{
									A(j, k) = aj[k];
									Q(j, k) = qj[k];
									A(i, k) = ai[k];
									Q(i, k) = qi[k];
								}
							}
						}
					}
					f[j] = tid + 1;
					// printf("tid:%d, f:%d, j:%d\n", tid, f[j], j);
				}
				f[M-1] = tid + 1;
				// printf("tid:%d, f:%d, j:%d\n", tid, f[M-1], M-1);
			}
		#pragma omp barrier
		}
	}

	// 串行
	if (THREADS_NUM == 1)
	{
		for (j = 0; j < M; j++)
			for (i = j + 1; i < M; i++)
			{
				sq = sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
				c = A(j, j) / sq;
				s = A(i, j) / sq;

				for (k = 0; k < M; k++)
				{
					aj[k] = c * A(j, k) + s * A(i, k);
					qj[k] = c * Q(j, k) + s * Q(i, k);
					ai[k] = (-s) * A(j, k) + c * A(i, k);
					qi[k] = (-s) * Q(j, k) + c * Q(i, k);
				}
				
				for (k = 0; k < M; k++)
				{
					A(j, k) = aj[k];
					Q(j, k) = qj[k];
					A(i, k) = ai[k];
					Q(i, k) = qi[k];
				}
			}
	}

	// 整理计算结果
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			R(i, j) = A(i, j);

	// 对Q进行转置
	for (i = 0; i < M; i++)
		for (j = i + 1; j < M; j++)
		{
			temp = Q(i, j);
			Q(i, j) = Q(j, i);
			Q(j, i) = temp;
		}

	end = omp_get_wtime();
	ParTime = (end - start) * 1e3;	//ms

	// 输出结果
	/*
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++) 
			printf("%f\t", A0(i, j));
		printf("\n");
	}*/

	//printf("\nOutput of QR operation\n");
	/*
	printf("\nMatrix R:\n");
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
			printf("%f\t", R(i, j));
		printf("\n");
	}*/
	
	/*
	printf("\nMatrix Q:\n");
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
			printf("%f\t", Q(i, j));
		printf("\n");
	}*/
	
	if (TestQR(A0, Q, R, M))
		printf("Success!\n");
	else
		printf("Singular Matrix!\n");
	printf("\nPar Time: %f\n", ParTime);

	return(0);
}
