#!/bin/bash

# Default number of processes
NUM_PROCESSES=${1:-4}
CITIES=${2:-10}
ITERATIONS=${3:-1000}

mpicc tspBruteForceParallelMPI.c -o tspBruteForceParallelMPI -lm && mpirun -np $NUM_PROCESSES ./tspBruteForceParallelMPI $CITIES $ITERATIONS