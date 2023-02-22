#!/bin/bash
echo "benchmark MPC"

SAMPLESDIR=../GenerateSamplesMATLAB/samples/rectangularSamples/
for j in {1..2}
do
	for elem in "FLOAT" "DOUBLE" 
	do
		for N in 2 3 4 6 8 16 32 64
		do
			#mkdir -p bins/${elem}/N${N}/
			for LS in 3
			do
				for i in {1..5}
				do
				echo
				echo $elem $N
				sample="${SAMPLESDIR}samplesMPC_N${N}.bin"
				./bins/${elem}/N${N}/benchmarkMPC_${elem}_LS${LS}_N${N} $sample bins/${elem}/N${N}/benchmarkMPC_${elem}_LS${LS}_N${N}.csv
				done
			done
		done
	done
done