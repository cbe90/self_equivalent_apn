To obtain the results for the 7-bit case, apply the following procedure:

1.) Compile the program with "g++ -O2 main_7bit.c -o 7bit" and run it.
For each of the 27 possible tuples, it will perform a quick check whether the tuple deserves further consideration or whether it can never yield an APN permutation. Set TIMEOUT to 900. Only the tuples corresponding to Classes 1, 16, 17, 18, 22, and 23 should remain. If more tuples remain, increase TIMEOUT.

2) Check those remaining 6 cases seperately by using the programs in the sub-directories.
