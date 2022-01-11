#ifndef UTIL_H
#define UTIL_H
#define N 1000


void init_array_d1(int A[], int n);
void init_array_d2(int A[][N]);
void campare_d1(const char *s, int A[], int B[], int n);
void campare_d2(const char *s, int A[][N], int B[][N]);

void init_array_d1_double(double A[], int n);
void init_array_d2_double(double A[][N]);
void campare_d1_double(const char *s, double A[], double B[], int n);
void campare_d2_double(const char *s, double A[][N], double B[][N]);

void print_array_d1(const char *s, int A[], int n);
#endif