#!/bin/bash

mkdir -p out
mkdir -p sol

# start all instances of the search for 7-bit shift-invariant APN permutations with sbox[1] = 3, sbox[3] = i
for i in {1..126}
do
	./det7_class1 sol/3_$i.txt 3 $i > out/3_$i.out &
done
