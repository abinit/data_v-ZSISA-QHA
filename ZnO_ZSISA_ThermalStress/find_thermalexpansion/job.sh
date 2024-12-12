#!/bin/bash

#SBATCH --job-name ZnO
#SBATCH --ntasks 40
#SBATCH --mem-per-cpu 2700mb
#SBATCH --output 132.log
#SBATCH --error 132.log
#SBATCH --time 0:20:00
##SBATCH --partition debug

export OMP_NUM_THREADS=1
export ABI_HOME=/home/ucl/modl/srostami/abinit-9.10.3 # Change this line
export PATH=$ABI_HOME/src/98_main/:$PATH      # Do not change this line: path to executable
export ABI_TESTS=$ABI_HOME/tests/            # Do not change this line: path to tests dir
export ABI_PSPDIR=/home/ucl/modl/srostami/Psps_for_tests/  # Do not change this line: path to pseudos dir

export ABI_SRC=$ABI_HOME/_build_/src/98_main/

module load releases/2021b intel libxc HDF5 netCDF netCDF-Fortran
mpirun $ABI_SRC/abinit Relax2.abi  > run.log 2> run.err
