#!/bin/bash

# Default number of cities
CITIES=${1:-10}
ITERATIONS=${2:-1000}

# Compile the tspBruteForce.c file
gcc tspBruteForce.c -o tspBruteForce -lm

# Run the compiled program with the specified number of cities
./tspBruteForce $CITIES $ITERATIONS