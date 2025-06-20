#!/bin/bash

# Default number of cities
CITIES=${1:-10}

# Compile the tspBruteForce.c file
gcc tspBruteForce.c -o tspBruteForce -lm

# Run the compiled program with the specified number of cities
./tspBruteForce $CITIES