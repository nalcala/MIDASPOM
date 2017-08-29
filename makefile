MPICC = mpicc
CC = gcc
CCFLAGS = -Wall -lcblas -latlas -lm -O3

all: MIDASPOMD_MPI MIDASPOMDs_MPI MIDASPOMDs2_MPI MIDASPOMf_MPI MIDASPOM_loss_MPI MIDASPOMD_dieoff_MPI MIDASPOM_future_MPI MIDASPOM MIDASPOM_dieoff MIDASPOM_loss MIDASPOM_future MIDASPOMD_impute_MPI 

MIDASPOMD_MPI: sources/main_MIDASPOMD_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMD_MPI.out sources/main_MIDASPOMD_MPI.c $(CCFLAGS)

MIDASPOMDs_MPI: sources/main_MIDASPOMDs_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMDs_MPI.out sources/main_MIDASPOMDs_MPI.c $(CCFLAGS)

MIDASPOMDs2_MPI: sources/main_MIDASPOMDs2_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMDs2_MPI.out sources/main_MIDASPOMDs2_MPI.c $(CCFLAGS)

MIDASPOMf_MPI: sources/main_MIDASPOMf_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMf_MPI.out sources/main_MIDASPOMf_MPI.c $(CCFLAGS)

MIDASPOMD_dieoff_MPI: sources/main_MIDASPOMD_dieoff_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMD_dieoff_MPI.out sources/main_MIDASPOMD_dieoff_MPI.c $(CCFLAGS)

MIDASPOM_loss_MPI: sources/main_MIDASPOM_loss_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_loss_MPI.out sources/main_MIDASPOM_loss_MPI.c $(CCFLAGS)

MIDASPOM_future_MPI: sources/main_MIDASPOM_future_MPI.c
	$(MPICC) -o bin_linux/MIDASPOM_future_MPI.out sources/main_MIDASPOM_future_MPI.c $(CCFLAGS)

MIDASPOMD_impute_MPI: sources/main_MIDASPOMD_impute_MPI.c
	$(MPICC) -o bin_linux/MIDASPOMD_impute_MPI.out sources/main_MIDASPOMD_impute_MPI.c $(CCFLAGS)

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
