This program finds all shift-invariant APN permutations in 7-bit up to affine equivalence.

To compile:
g++ -O2 main_det7_class1.c -o det7_class1

To complete the search, we need to start 5*126 = 630 processes. The following scripts starts one chunk of 126 processes each (in background):
./start_7_class1_1.sh
./start_7_class1_2.sh
./start_7_class1_3.sh
./start_7_class1_4.sh
./start_7_class1_5.sh

To kill all the processes, run:
./kill_7_class1.sh


The .txt files will contain all the shift-invariant APN permutations, though the same solution (up to equivalence) can be found multiple times.
The .out files show the console output of each process.
