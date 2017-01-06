# MIDASPOM
Metapopulation Inference from Data in Stochastic Patch Occupancy Models (MIDASPOM) is a Bayesian Inference method for patch occupancy (presence/absence) data.

## Install and run
On a Linux system, put the downloaded file into a suitable directory.   

Next, unzip the file. Type the command 
```bash
unzip MIDASPOM.zip
```
in a terminal or right click and select "Extract" in the menu. 

Move into the ‘MIDASPOM’ directory 
```bash
cd MIDASPOM
```
You can run the software by typing
```bash
./bin/MIDASPOM.out
```
in the ‘MIDASPOM’ directory, and adding the appropriate flags (e.g., "-i occupancies.txt" to use file occupancies.txt as input data). The file "Manual.pdf" from directory "manual" provides a complete description of software usage with a list of all possible flags. There is also a parallel processing versio of the program that you can run by typing 
```bash
./bin/MIDASPOM_MPI.out
```

File "run_examples.sh" is a script that runs the software on a few examples, and produces plots using the R scripts from folder "Rscript". Just type 
```bash
./run_examples.sh
```
in a terminal to compute these examples. The output files and plots are written in the "examples/output" folder. Type 
```bash
./run_examples_MPI.sh
```
to run the MPI version of the program for faster execution time.

Notes:
1. In order to compile the sources yourself, you will need the gcc compiler, the CBLAS librairies (from the ATLAS package), and the mpicc compiler and the MPI librairies to compile the parallel versions of the sources. To compile, you can type "make" in a terminal in the MIDASPOM directory  
2. In order to launch the script compute.R, you will need the software R, available from the CRAN website https://cran.r-project.org/

