#!/bin/bash
#PBS -N zTraclp
#PBS -A UWAS0108
#PBS -l walltime=00:15:00
#PBS -q develop
#PBS -j oe
#PBS -k eod
#PBS -m be
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=8

module load cdo
 
run=aa_lpfrz
run_dir=/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/dpscream_rce_large_3km_${run}/run/
ds_file=$run_dir/dpscream_rce_large_3km_${run}.eam.h0.2000-01-01-00000.nc

cdo -seltimestep,-20/-1 $ds_file $run_dir/dpscream_rce_large_3km_${run}_h0_last5days.nc

# cdo -fldmean -timmean -seltimestep,-20/-1 -selvar,Z3 $ds_file $run_dir/dpscream_rce_large_3km_${run}_Z3_mean.nc

# tracer=BCU
# tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
# cdo -mulc,-1 -ln -selvar,$tracer $ds_file $tracer_file # -seltimestep,-20/-1 

# tracer=NUC
# tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
# cdo -mulc,-1 -ln -selvar,$tracer $ds_file $tracer_file # -seltimestep,-20/-1

# tracer=W_NUC
# tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
# cdo -setname,$tracer -div -selvar,$tracer $ds_file -selvar,NUC $ds_file $tracer_file

# tracer=NI_NUC
# tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
# cdo -setname,$tracer -div -selvar,$tracer $ds_file -selvar,NUC $ds_file $tracer_file


echo "done"
