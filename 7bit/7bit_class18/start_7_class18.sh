#!/bin/bash

mkdir -p out
mkdir -p sol

# start all instances of the search for 7-bit APN permutations in class 18 with sbox[1] = 1 or sbox[1] = 8
./det7_class18 sol/1.txt 1 > out/1.out &
./det7_class18 sol/8.txt 8 > out/8.out &
