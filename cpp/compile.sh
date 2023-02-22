#!/bin/bash
echo "compile MPC"

make clean

for elem in "FLOAT" "DOUBLE" 
do
	for N in 2 3 4 6 8 16 32 64
	do
		mkdir -p bins/${elem}/N${N}/
		for LS in 3
		do
			make DEFINES="-DELEM_${elem} -DN_HOR=${N}  -DLS=${LS}"
			mv build/bin/runner bins/${elem}/N${N}/benchmarkMPC_${elem}_LS${LS}_N${N}
			make clean
		done
	done
done