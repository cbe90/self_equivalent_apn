To obtain the results for the 8-bit case, apply the following procedure:

1.) Compile the program with "g++ -O2 main_8bit.c -o 8bit" and run it.
For each of the 32 possible tuples, it will perform a quick check whether the tuple deserves further consideration or whether it can never yield an APN permutation. Set TIMEOUT to 900. Only the tuples corresponding to Classes 1,2,22, and 30 should remain. If more tuples remain, increase TIMEOUT.

2) Check those remaining four cases seperately by using the programs in the sub-directories.
