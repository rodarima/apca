#!/bin/bash

if [ "$1" == "" ]; then
	echo "Usage: $0 <seconds>"
	echo ""
	echo "Dummy work during the especified time."
	exit 1
fi

#pwd
#cat /proc/cpuinfo
#uname -a
#ip addr
#hostname
yes > /dev/null &
PID=$!
sleep $1
kill -9 $PID
echo "$SLURM_PROCID:OK"
