#!/bin/bash

mkdir -p out
mkdir -p sol

# start 32 instances of the search and run them for 7 days
for i in {1..32}
do
	./rand7_class16 sol/$i.txt 604800 > out/$i.out &
done
