#!/bin/bash

# Default number of processes
NUM_PROCESSES=${1:-4}
CITIES=${2:-10}
ITERATIONS=${3:-1000}
NUM_THREADS=${4:-8}

mpicc -fopenmp tspBruteForceParallelHybrid.c -o tspBruteForceParallelHybrid -lm

mpirun -np $NUM_PROCESSES ./tspBruteForceParallelHybrid $CITIES $ITERATIONS $NUM_THREADS