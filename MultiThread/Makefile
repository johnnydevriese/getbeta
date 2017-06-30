CC = g++
FC = gfortran
OBJS = myompsubs.o mychebfuns.o HamOpsubs.o ChebPolysubs.o ChebPolyOpsubs.o myfortompsubs.o
LIBS = -llapack -lblas -lm -lfftw3 -larpack -lfftw3_omp -lgfortran -lgsl -lgslcblas
DIR = -I/home/jared/Dropbox/MultiThread/include -I/usr/local/arpack++/include
FLAGS = -O3 -fopenmp

all: test

test: ${OBJS} test.cpp
	${CC} ${FLAGS} ${DIR} -o test test.cpp ${OBJS} ${LIBS}

myompsubs.o: myompsubs.cpp
	${CC} ${FLAGS} ${DIR} -c myompsubs.cpp ${LIBS}

mychebfuns.o: mychebfuns.cpp
	${CC} ${FLAGS} ${DIR} -c mychebfuns.cpp ${LIBS}

HamOpsubs.o: HamOpsubs.cpp 
	${CC} ${FLAGS} ${DIR} -c HamOpsubs.cpp ${LIBS}

ChebPolysubs.o: ChebPolysubs.cpp
	${CC} ${FLAGS} ${DIR} -c ChebPolysubs.cpp ${LIBS}

ChebPolyOpsubs.o: ChebPolyOpsubs.cpp
	${CC} ${FLAGS} ${DIR} -c ChebPolyOpsubs.cpp ${LIBS}

myfortompsubs.o: myfortompsubs.f90
	${FC} ${FLAGS} ${DIR} -c myfortompsubs.f90 ${LIBS}
 
clean:
	rm test *.o
