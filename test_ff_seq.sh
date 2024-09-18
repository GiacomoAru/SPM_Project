#!/bin/bash
if [ $# -lt 2 ]; then
    echo "Usage: $0 MACHINE_ID TOT_ITERATION"
    exit 1
fi

MACHINE_ID=$1
NUM_CORES=$2

THREADS_STEP=2
REPETITIONS=5
DEFAULT_N=1250
STEP_N=64
STEP_2N=64000

LOG_FILE_DM=./log/seq_dm.csv
LOG_FILE_SM=./log/seq_sm.csv

PROG_DM=./sequential_dm
PROG_SM=./sequential_sm

######################## TESTING STRONG SCALABILITY ############################

echo "Computing sequential time reference"
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] N = $DEFAULT_N, m = $MACHINE_ID, l = $LOG_FILE_SM"
    $PROG_SM -N $DEFAULT_N -m $MACHINE_ID -l $LOG_FILE_SM
done
for rep in $(seq 1 $REPETITIONS); do
    echo "[$rep/$REPETITIONS] N = $DEFAULT_N, m = $MACHINE_ID, l = $LOG_FILE_DM"
    $PROG_DM -N $DEFAULT_N -m $MACHINE_ID -l $LOG_FILE_DM
done

######################## QUADRATIC GROWTH ##############################

echo "Weak scalability comparison - single matrix - quadratic growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    step=$((t == 0 ? t + 1 : t))
    n=$((step*STEP_N))
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] N = $n, m = $MACHINE_ID, l = $LOG_FILE_SM"
        $PROG_SM -N $n -m $MACHINE_ID -l $LOG_FILE_SM
    done
done

echo "Weak scalability comparison - double matrix - quadratic growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    step=$((t == 0 ? t + 1 : t))
    n=$((step*STEP_N))
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] N = $n, m = $MACHINE_ID, l = $LOG_FILE_DM"
        $PROG_DM -N $n -m $MACHINE_ID -l $LOG_FILE_DM
    done
done

######################## LINEAR GROWTH  ##############################

echo "Weak scalability comparison - single matrix - linear growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    step=$((t == 0 ? t + 1 : t))
    n=$((step*STEP_2N))
    sq_n=$(awk "BEGIN {print sqrt($n)}")
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] N = $sq_n, m = $MACHINE_ID, l = $LOG_FILE_SM"
        $PROG_SM -N $sq_n -m $MACHINE_ID -l $LOG_FILE_SM
    done
done

echo "Weak scalability comparison - double matrix - linear growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    step=$((t == 0 ? t + 1 : t))
    n=$((step*STEP_2N))
    sq_n=$(awk "BEGIN {print sqrt($n)}")
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] N = $sq_n, m = $MACHINE_ID, l = $LOG_FILE_DM"
        $PROG_DM -N $sq_n -m $MACHINE_ID -l $LOG_FILE_DM
    done
done

