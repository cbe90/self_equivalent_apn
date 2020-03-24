This program searches for APN permutations in 7-bit which admit an automorphism from class 18, up to affine equivalence.

To compile:
g++ -O2 main_det7_class18.c -o det7_class18

To complete the search, we need to start 2 processes. The following scripts start them (in background):
./start_7_class18.sh

To kill all the processes, run:
./kill_7_class18.sh


The .txt files will contain all the APN permutations in class 18, though the same solution (up to equivalence) can be found multiple times.
The .out files show the console output of each process.
