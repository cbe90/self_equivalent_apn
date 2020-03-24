#!/bin/bash

mkdir -p out
mkdir -p sol

# start all instances of the search for 7-bit shift-invariant APN permutations with sbox[1] = 0xb, sbox[3] = i
for i in {1..126}
do
	./det7_class1 sol/11_$i.txt 11 $i > out/11_$i.out &
done
