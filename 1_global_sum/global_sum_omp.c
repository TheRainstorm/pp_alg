#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <omp.h>

int op(int a, int b){
    long long loop = 1000000000;
    while(loop--);
    // sleep(1);
    return a + b;
}

void basic_sum(int *data, int n){
    int res = 0;
    for(int i=0; i<n; i++){
        res = op(res, data[i]);
    }
    
    for(int i=0; i<n; i++){
        data[i] = res;
    }
}

void BiTreeSum(int *data, int n) {
    //Thanks to original author: SA21严雪文
	int i, step;

	for (step = 1; step != n; step = step * 2) {	// 每轮计算的节点是step的倍数
		#pragma omp parallel for shared(data) private(i)
		for (i = 0; i < n; i++) {
			if(!(i % (step * 2))) 	// 下一轮中留下来的节点计算结果
            {
                // int tid = omp_get_thread_num();
                // printf("step: %d tid: %d i: %d, add: %d\n",step, tid, i, i + step);
				data[i] = op(data[i], data[i + step]);
            }
		}
	}

    for(int i=1; i<n; i++){
        data[i] = data[0];
    }
}

void BiTreeSum2(int *data, int n) {
	for (int step = 1; step != n; step = step * 2) {	// 每轮计算的节点是step的倍数
		#pragma omp parallel shared(data)
        {
            int i = omp_get_thread_num();
            if(!(i % (step * 2))) 	// 下一轮中留下来的节点计算结果
            {
                // printf("step: %d tid: %d i: %d, add: %d\n",step, tid, i, i + step);
                data[i] = op(data[i], data[i + step]);
            }
        }
	}

    for(int i=1; i<n; i++){
        data[i] = data[0];
    }
}

void ButterflySum(int *data, int n) {
    //Thanks to original author: SA21严雪文
	int i, step = 1;
    int *temp = (int *)malloc(n*sizeof(int));
	while (step < n) {
		#pragma omp parallel for shared(data, temp) private(i)
		for (i = 0; i < n; i++) {
            temp[i] = op(data[i], data[i ^ step]);
		}
        
		#pragma omp parallel for shared(data, temp) private(i)
		for (i = 0; i < n; i++)
			data[i] = temp[i];

		step <<= 1;
	}
}

void ButterflySum2(int *data, int n) {
	int i, step = 1;
    int *temp = (int *)malloc(n*sizeof(int));
	while (step < n) {
		#pragma omp parallel shared(data, temp)
        {
            int i = omp_get_thread_num();
            // printf("step: %d tid: %d i: %d, add: %d\n",step, tid, i, i ^ step);
            temp[i] = op(data[i], data[i ^ step]);

            // #pragma omp barrier
            data[i] = temp[i];
        }

		step <<= 1;
	}
}

int main(int argc, char *argv[]){
    int select = 0;
    switch (argc){
    case 2:
        select = atoi(argv[1]);
        break;
    case 1:
        break;
    default:
        printf("Usage: program 0|1\n");
        exit(0);
    }

    int nthread = omp_get_max_threads();
    printf("nthread: %d\n", nthread);
    
    int *data = (int *)malloc(nthread*sizeof(int));
	#pragma omp parallel
	{
		int id;
		id = omp_get_thread_num();
		data[id] = id;
	}

    long t = -clock();
    switch (select){
    case 4:
        BiTreeSum2(data, nthread);
        break;
    case 3:
        ButterflySum2(data, nthread);
        break;
    case 2:
        BiTreeSum(data, nthread);
        break;
    case 1:
        ButterflySum(data, nthread);
        break;
    case 0:
        basic_sum(data, nthread);
    default:
        break;
    }
    t += clock();

    printf("sum \n");
    for (int i = 0; i < nthread; i++)
		printf("%d ", data[i]);
	printf("\n");
    
    printf("time elapsed: %f\n", (double)t/CLOCKS_PER_SEC);
    return 0;
}

/*
output:

➜  1_global_sum time ./global_sum_omp 1
nthread: 8
sum
28 28 28 28 28 28 28 28
time elapsed: 101.743557
./global_sum_omp 1  101.73s user 0.02s system 747% cpu 13.606 total
➜  1_global_sum time ./global_sum_omp 0
nthread: 8
sum
28 28 28 28 28 28 28 28
time elapsed: 21.606844
./global_sum_omp 0  21.60s user 0.01s system 100% cpu 21.567 total


21.567 / 13.606 = 1.585的加速比
理论加速比：basic_sum需要进行7次op操作，而ButterflySum平均每个线程需要进行3次op操作，
因此，理论加速比为7/3 =2.33

Queston：有什么profile工具，可以看到ButterflySum两个for循环的执行时间呢？
*/