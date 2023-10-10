#!/bin/bash

# Compile 3d_matrix.cpp
g++ -I /usr/include/eigen3 -O3 3d_matrix.cpp -o matrix

# Execute matrix
for ((n = 3; n <= 30; n++))
do
    ./matrix $n
done

# Plot the results
python3 plot_mem.py
