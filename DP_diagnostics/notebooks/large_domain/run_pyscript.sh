#!/bin/bash -l
#PBS -N py_mhist_plot
#PBS -A UWAS0108
#PBS -j oe
#PBS -k eod
#PBS -q casper
#PBS -M smturbev@uw.edu
#PBS -m be
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mem=10GB

export TMPDIR=${SCRATCH}/temp
mkdir -p ${TMPDIR}

### Load Conda/Python module and activate NPL environment

module load conda
conda activate npl-2023a

### Run analysis script
python plot_mhist.py

echo "done"
