CC=mpicc  # variable donnant le compilateur
CFLAGS=-g --std=c99 -pedantic -Wall -Wextra -Wmissing-prototypes

all: project2_DONVAL_CHIODO_main.o project2_DONVAL_CHIODO_map.o project2_DONVAL_CHIODO_explicit.o project2_DONVAL_CHIODO_implicit.o
	$(CC) project2_DONVAL_CHIODO_main.o project2_DONVAL_CHIODO_map.o project2_DONVAL_CHIODO_explicit.o project2_DONVAL_CHIODO_implicit.o -lm -o projet2

project2_DONVAL_CHIODO_main.o: project2_DONVAL_CHIODO_main.c project2_DONVAL_CHIODO_main.h
	$(CC) $(CFLAGS) -c project2_DONVAL_CHIODO_main.c 

project2_DONVAL_CHIODO_map.o: project2_DONVAL_CHIODO_map.c project2_DONVAL_CHIODO_map.h
	$(CC) $(CFLAGS) -c project2_DONVAL_CHIODO_map.c 

project2_DONVAL_CHIODO_explicit.o: project2_DONVAL_CHIODO_explicit.c project2_DONVAL_CHIODO_explicit.h
	$(CC) $(CFLAGS) -c project2_DONVAL_CHIODO_explicit.c 

project2_DONVAL_CHIODO_implicit.o: project2_DONVAL_CHIODO_implicit.c project2_DONVAL_CHIODO_implicit.h
	$(CC) $(CFLAGS) -c project2_DONVAL_CHIODO_implicit.c 

runSriLanka:
	./projet2 param.txt hw2/sriLanka.dat 0

runReflection:
	mpirun -np 1 ./projet2 param.txt hw2/reflection.dat 0
clean:
	rm *.o projet2
