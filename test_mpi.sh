#!/bin/sh
#SBATCH -p normal
#SBATCH -N 8
#SBATCH --ntasks-per-node=1
#SBATCH -t 01:00:00

MACHINE_ID=3
NUM_PROC=8

REPETITIONS=3
DEFAULT_N=2048
STEP_N=512

LOG_FILE=./log/mpi.csv

PROG=./mpi_wavefront


######################## TESTING STRONG SCALABILITY ############################
echo "Testing strong scalability"
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] P = 1, W = 2, N = $DEFAULT_N, G = 256, m = $MACHINE_ID, l = $LOG_FILE"
    mpirun -n 1 $PROG -W 2 -N $DEFAULT_N -G 256 -m $MACHINE_ID -l $LOG_FILE
done
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] P = 2, W = 2, N = $DEFAULT_N, G = 64, m = $MACHINE_ID, l = $LOG_FILE"
    mpirun -n 2 $PROG -W 2 -N $DEFAULT_N -G 64 -m $MACHINE_ID -l $LOG_FILE
done
for t in $(seq 3 1 $NUM_PROC); do
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] P = $t, W = 2, N = $DEFAULT_N, G = 32, m = $MACHINE_ID, l = $LOG_FILE"
        mpirun -n $t $PROG -W 2 -N $DEFAULT_N -G 32 -m $MACHINE_ID -l $LOG_FILE
    done
done


######################## QUADRATIC GROWTH ##############################
echo "Testing weak scalability - quadratic growth"
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] P = 1, W = 2, N = 512, G = 256, m = $MACHINE_ID, l = $LOG_FILE"
    mpirun -n 1 $PROG -W 2 -N 512 -G 256 -m $MACHINE_ID -l $LOG_FILE
done
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] P = 2, W = 2, N = 1024, G = 64, m = $MACHINE_ID, l = $LOG_FILE"
    mpirun -n 2 $PROG -W 2 -N 1024 -G 64 -m $MACHINE_ID -l $LOG_FILE
done
for t in $(seq 3 1 $NUM_PROC); do
    n=$((t*STEP_N))
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] P = $t, W = 2, N = $n, G = 32, m = $MACHINE_ID, l = $LOG_FILE"
        mpirun -n $t $PROG -W 2 -N $n -G 32 -m $MACHINE_ID -l $LOG_FILE
    done
done