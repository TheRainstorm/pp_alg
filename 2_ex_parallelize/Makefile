
LIBS = -L. -lutil
CFLAGS = -O2 -fopenmp $(LIBS)
CFLAGS2 = -O2 $(LIBS)
OMP_FLGA = -fopenmp 
VEC_FLAG = -ftree-loop-vectorize -march=core-avx2 -fopt-info-vec

EXE = ex_1_3 ex_1_3_2 ex_3_2 ex_3_3 ex_3_5_1 ex_3_5_2 ex_3_5_3 test_F
all: $(EXE)

ex_1_3: libutil.a ex_1_3.c
	gcc ex_1_3.c $(CFLAGS) -o $@
	gcc ex_1_3.c $(CFLAGS) $(VEC_FLAG) -o ex_1_3_vec

ex_1_3_2: libutil.a ex_1_3_2.c
	gcc ex_1_3_2.c $(CFLAGS) -o $@
	# gcc ex_1_3_2.c $(CFLAGS) $(VEC_FLAG) -o ex_1_3_2_vec

ex_3_2: libutil.a ex_3_2.c
	gcc ex_3_2.c -O2 -L. -lutil -o ex_3_2_nomp
	gcc ex_3_2.c -O2 -L. -lutil $(VEC_FLAG) -o ex_3_2_nomp_vec
	gcc ex_3_2.c $(CFLAGS) -o $@
	gcc ex_3_2.c $(CFLAGS) $(VEC_FLAG) -o ex_3_2_vec

ex_3_3: libutil.a ex_3_3.c
	gcc ex_3_3.c $(CFLAGS) -o $@
	# gcc ex_3_3.c $(CFLAGS) $(VEC_FLAG) -o ex_3_3_vec

ex_3_5_1: libutil.a ex_3_5_1.c
	gcc ex_3_5_1.c $(CFLAGS2) -o $@
	gcc ex_3_5_1.c $(CFLAGS2) $(OMP_FLGA) -o ex_3_5_1_omp
	gcc ex_3_5_1.c $(CFLAGS2) $(OMP_FLGA) $(VEC_FLAG) -o ex_3_5_1_vec

ex_3_5_2: libutil.a ex_3_5_2.c
	gcc ex_3_5_2.c $(CFLAGS2) -g -o $@
	gcc ex_3_5_2.c $(CFLAGS2) $(OMP_FLGA) -o ex_3_5_2_omp
	gcc ex_3_5_2.c $(CFLAGS2) $(OMP_FLGA) $(VEC_FLAG) -o ex_3_5_2_vec

ex_3_5_3: libutil.a ex_3_5_3.c
	gcc ex_3_5_3.c $(CFLAGS2) -g -o $@
	gcc ex_3_5_3.c $(CFLAGS2) $(OMP_FLGA) -o ex_3_5_3_omp
	gcc ex_3_5_3.c $(CFLAGS2) $(OMP_FLGA) $(VEC_FLAG) -o ex_3_5_3_vec

test_F: libutil.a test_F.c
	gcc test_F.c $(CFLAGS) -g -o $@

libutil.a: util.c util.h
	gcc -c util.c
	ar rcs $@ util.o

clean:
	rm $(EXE)
	rm *.o
	rm *.a
	rm *_omp
	rm *_nomp
	rm *_vec