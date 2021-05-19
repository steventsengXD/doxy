# Include this so that compile/run errors terminate script
set -e

# Place all necessary fortran programs/modules into a single variable
myprogramfiles="write_netcdf.f90 command_line.f90 random_mod.f90 sweetener.f90 ising.f90"

# Name of the compiled file
outfile="ising"

# Name of compiler
fc=gfortran

# We use nf-config to grab the compile and link flags, backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`

# Command to compile program/modules along with necessary flags
$fc -g -std=f2008 -Wall $fflags $myprogramfiles $flibs -o $outfile

# Terminal command to run the Ising model with an initial random state, grid of size
# 50x50, over 1000000 timesteps, temperature = 2 and interaction strength 1.
# These are the same values that the program will default to if there are no/incorrect
# command line inputs.
./ising N=50 T=1000000 init=F beta=0.5 J=1
# ./ising

# Plotting visualizations of data in python3
python3 ising.py
