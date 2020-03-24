To obtain the results for the 6-bit case, apply the following procedure:

1.) Compile the program with "g++ -O2 main_6bit.c -o 6bit" and run it.
For each of the 17 possible tuples, it will perform a quick check whether the tuple deserves further consideration or whether it can never yield an APN permutation. Set TIMEOUT to 900. Only the tuples corresponding to Class 5 (A = B = Comp(X^6+X^5+X^4+X^3+X^2+X+1)) and Class 14 (A = B = Comp(X^2+1) \oplus Comp(X^2+1) \oplus Comp(X^2+1)) should remain. If more tuples remain, increase TIMEOUT.

2) Check those remaining two cases seperately by using the programs in the "/6bit_class5" and "/6bit_class14" directories.
