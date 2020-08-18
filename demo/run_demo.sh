#!/bin/bash

#------------------------------------------------------------------------------|
#   Run fortran demo of hot-tail calculation
#------------------------------------------------------------------------------|
echo "Running fortran demo..."
./hot_tail_demo
mv hot_tails.dat dat/hot_tails_fortran.dat

#------------------------------------------------------------------------------|
#   Repeat calculation for same plasma parameters with python
#------------------------------------------------------------------------------|
echo "Running python demo..."
python scripts/hot_tails.py

#------------------------------------------------------------------------------|
#   Repeat calculation for same plasma parameters with MATLAB
#------------------------------------------------------------------------------|
echo "Running MATLAB demo..."
matlab -nodesktop -nosplash -noFigureWindow -r \
    'try, run ("scripts/hot_tails.m"); end; quit;' > /dev/null

#------------------------------------------------------------------------------|
#   Run comparison .py-script
#------------------------------------------------------------------------------|
echo "Compare results from each calculation..."
python scripts/compare_populations.py

#------------------------------------------------------------------------------|
