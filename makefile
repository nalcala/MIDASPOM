MPICC = mpicc
CC = gcc
CCFLAGS = -latlas -lcblas -lm -O3

all: MIDASPOM_MPI MIDASPOM_loss_MPI MIDASPOM_dieoff_MPI MIDASPOM_future_MPI MIDASPOM MIDASPOM_dieoff MIDASPOM_loss MIDASPOM_future

MIDASPOM_MPI: sources/main_MIDASPOM_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_MPI.out sources/main_MIDASPOM_MPI.c $(CCFLAGS)

MIDASPOM_dieoff_MPI: sources/main_MIDASPOM_dieoff_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_dieoff_MPI.out sources/main_MIDASPOM_dieoff_MPI.c $(CCFLAGS)

MIDASPOM_loss_MPI: sources/main_MIDASPOM_loss_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_loss_MPI.out sources/main_MIDASPOM_loss_MPI.c $(CCFLAGS)

MIDASPOM_future_MPI: sources/main_MIDASPOM_future_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_future_MPI.out sources/main_MIDASPOM_future_MPI.c $(CCFLAGS)

MIDASPOM: sources/main_MIDASPOM.c
	$(CC) -o bin_linux/MIDASPOM.out sources/main_MIDASPOM.c $(CCFLAGS)

MIDASPOM_dieoff: sources/main_MIDASPOM_dieoff.c
	$(CC) -o bin_linux/MIDASPOM_dieoff.out sources/main_MIDASPOM_dieoff.c $(CCFLAGS)

MIDASPOM_loss: sources/main_MIDASPOM_loss.c
	$(CC) -o bin_linux/MIDASPOM_loss.out sources/main_MIDASPOM_loss.c $(CCFLAGS)

MIDASPOM_future: sources/main_MIDASPOM_future.c
	$(CC) -o bin_linux/MIDASPOM_future.out sources/main_MIDASPOM_future.c $(CCFLAGS)

clean:
	rm sources/*.o
