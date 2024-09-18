#!/bin/bash
if [ $# -lt 3 ]; then
    echo "Usage: $0 MACHINE_ID TOT_ITERATION MAX_GRANULARITY"
    exit 1
fi

MACHINE_ID=$1
NUM_CORES=$2
MAX_GRANULARITY=$3

THREADS_STEP=2
GRAIN_BASE=2
REPETITIONS=5
DEFAULT_N=2000

LOG_FILE_DM=./log/ff_s_dm_g.csv
LOG_FILE_SM=./log/ff_s_sm_g.csv

PROG_DM=./fastflow_stream_dm
PROG_SM=./fastflow_stream_sm

######################## TESTING STRONG SCALABILITY ############################
echo "Testing granularity - single matrix"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    for g in $(seq 0 1 $MAX_GRANULARITY); do
        gran=$(awk "BEGIN { print exp($g * log(2)) }")
        for rep in $(seq 1 $REPETITIONS); do
            echo "[$rep/$REPETITIONS] W = $thr, N = $DEFAULT_N, G = $gran, m = $MACHINE_ID, l = $LOG_FILE_SM"
            $PROG_SM -W $thr -N $DEFAULT_N -G $gran -m $MACHINE_ID -l $LOG_FILE_SM
        done
    done
done
echo "Testing granularity - double matrix"
for t in $(seq 0 $THREADS_STEP $NUM_CORES); do
    thr=$((t == 0 ? t + 1 : t))
    for g in $(seq 0 1 $MAX_GRANULARITY); do
        gran=$(awk "BEGIN { print exp($g * log(2)) }")
        for rep in $(seq 1 $REPETITIONS); do
            echo "[$rep/$REPETITIONS] W = $thr, N = $DEFAULT_N, G = $gran, m = $MACHINE_ID, l = $LOG_FILE_DM"
            $PROG_DM -W $thr -N $DEFAULT_N -G $gran -m $MACHINE_ID -l $LOG_FILE_DM
        done
    done
done
