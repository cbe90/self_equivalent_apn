This program finds all APN permutations in 8-bit from class 1 up to affine equivalence.

To compile:
g++ -O2 main_det8_class1.c -o det8_class1

To complete the search, we need to start 255 processes. The following script starts them all (in background):
./start_8_class1.sh

To kill all the processes, run:
./kill_8_class1.sh


The .txt files will contain all the 8-bit APN permutations in class 1, though the same solution (up to equivalence) can be found multiple times.
The .out files show the console output of each process.
