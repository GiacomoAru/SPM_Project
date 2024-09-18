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

LOG_FILE_DM=./log/ff_d_dm.csv
LOG_FILE_SM=./log/ff_d_sm.csv

PROG_DM=./fastflow_data_dm
PROG_SM=./fastflow_data_sm

######################## TESTING STRONG SCALABILITY ############################
echo "Testing strong scalability - single matrix"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    for rep in $(seq 1 $REPETITIONS); do  
        echo "[$rep/$REPETITIONS] W = $thr, N = $DEFAULT_N, m = $MACHINE_ID, l = $LOG_FILE_SM"
        $PROG_SM -W $thr -N $DEFAULT_N -m $MACHINE_ID -l $LOG_FILE_SM
    done
done
echo "Testing strong scalability - double matrix"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    for rep in $(seq 1 $REPETITIONS); do  
        echo "[$rep/$REPETITIONS] W = $thr, N = $DEFAULT_N, m = $MACHINE_ID, l = $LOG_FILE_DM"
        $PROG_DM -W $thr -N $DEFAULT_N -m $MACHINE_ID -l $LOG_FILE_DM
    done
done
######################## QUADRATIC GROWTH ##############################

echo "Testing weak scalability - single matrix - quadratic growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    n=$((thr*STEP_N))
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] W = $thr, N = $n, m = $MACHINE_ID, l = $LOG_FILE_SM"
        $PROG_SM -W $thr -N $n -m $MACHINE_ID -l $LOG_FILE_SM
    done
done

echo "Testing weak scalability - double matrix - quadratic growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    n=$((thr*STEP_N))
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] W = $thr, N = $n, m = $MACHINE_ID, l = $LOG_FILE_DM"
        $PROG_DM -W $thr -N $n -m $MACHINE_ID -l $LOG_FILE_DM
    done
done

######################## LINEAR GROWTH  ##############################

echo "Testing weak scalability - single matrix - linear growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    n=$((thr*STEP_2N))
    sq_n=$(awk "BEGIN {print sqrt($n)}")
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] W = $thr, N = $sq_n, m = $MACHINE_ID, l = $LOG_FILE_SM"
        $PROG_DM -W $thr -N $sq_n -m $MACHINE_ID -l $LOG_FILE_SM
    done
done

echo "Testing weak scalability - double matrix - linear growth"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    n=$((thr*STEP_2N))
    sq_n=$(awk "BEGIN {print sqrt($n)}")
    for rep in $(seq 1 $REPETITIONS); do
        echo "[$rep/$REPETITIONS] W = $thr, N = $sq_n, m = $MACHINE_ID, l = $LOG_FILE_DM"
        $PROG_DM -W $thr -N $sq_n -m $MACHINE_ID -l $LOG_FILE_DM
    done
done

