#!/bin/bash

mkdir -p out
mkdir -p sol

# start all instances of the search for 8-bit APN permutations in class 1 with sbox[1] = 1, sbox[3] = i
for i in {1..255}
do
	./det8_class1 sol/1_$i.txt 1 $i > out/1_$i.out &
done
