#!/bin/bash

# Se viene fornito un argomento, usa quello, altrimenti usa il valore di default
matrix_size="${1:-10}"
worker_num="${2:-2}"
granularity="${3:-2}"

echo "TEST with matrix_size = $matrix_size worker_num = $worker_num granularity = $granularity"

echo "TEST sequential_sm"
./sequential_dm -N $matrix_size

echo "TEST sequential_dm"
./sequential_sm -N $matrix_size

echo "TEST fastflow_data_sm"
./fastflow_data_sm -N $matrix_size -W $worker_num

echo "TEST fastflow_data_dm"
./fastflow_data_dm -N $matrix_size -W $worker_num

echo "TEST fastflow_stream_sm"
./fastflow_stream_sm -N $matrix_size -W $worker_num -G $granularity

echo "TEST fastflow_stream_dm"
./fastflow_stream_dm -N $matrix_size -W $worker_num -G $granularity

echo "TEST mpi_wavefront not mprun"
./mpi_wavefront -W 1 -N $matrix_size -G $granularity

echo "TEST mpi_wavefront"
mpirun -n $worker_num ./mpi_wavefront -W 1 -N $matrix_size -G $granularity