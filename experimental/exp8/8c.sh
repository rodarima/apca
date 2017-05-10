#!/bin/bash

RUNS=100000
#STEPS=$(echo $(seq 1 9; seq 10 10 100) | tr '\n' ' ')
STEPS=$(seq -s ' ' 1 9 | tr '\n' ' '; seq -s' ' 10 10 100)
BITS=10

STEPSFILE="data/8c-steps.csv"

test -f "$STEPSFILE" || echo "$STEPS" | tr ' ' '\n' > "$STEPSFILE"

for S in $STEPS; do

	OUTFILE="data/8c-$RUNS-$S-$BITS.csv"
	if [ -f "$OUTFILE" ]; then
		echo "Skipping file $OUTFILE. Already computed."
		continue
	fi
	echo "Running: ./8c $RUNS $S $BITS"
	echo "Computing file $OUTFILE"
	#./8c $RUNS $S $BITS > "$OUTFILE"

done
