#!/bin/bash


for kAngle in 200 300 400 500
do
	mkdir ./md/periodic-nonperiodic-comparison4/kAngle=${kAngle}
	for mol in 0 1 2 3 4 5 6 7 8 9
	do
		./main -config ./config/my.config -set main.lammps-kAngle.start ${kAngle} -set io.lammps-output ./md/periodic-nonperiodic-comparison4/kAngle=${kAngle}/ -set io.input ./data/periodic-nonperiodic-comparison/run_N=1054_nPos=86_nNeg=82_layers=5_molecules/molecule${mol} -set main.inputFile molecule${mol}.sim -set main.logFile molecule${mol}-log.sim -set main.dumpFile molecule${mol}.xyz -set main.topologyFile molecule${mol}.0time
	done
done

