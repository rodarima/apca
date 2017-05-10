#!/bin/bash

RUNS=100000
STEPS=$(echo $(seq 1 9; seq 10 10 100) | tr '\n' ' ')
BITS=10

for S in $STEPS; do

	echo "Running: ./exp8c $RUNS $S $BITS"
	./exp8c $RUNS $S $BITS > ../data/exp8c-$RUNS-$S-$BITS.csv

done
