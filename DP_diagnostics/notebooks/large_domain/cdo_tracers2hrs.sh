#!/bin/bash -v
#PBS -N zTracBranch
#PBS -A UWAS0108
#PBS -l walltime=00:30:00
#PBS -q develop
#PBS -j oe
#PBS -k eod
#PBS -m be
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=87

module load cdo
 
run=lpls300k_wbranch
run_dir=/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/dpscream_rce_large_3km_${run}/run/
ds_file=$run_dir/dpscream_rce_large_3km_${run}.eam.h0.2000-02-20-10800.nc
ds5day_file=$ds_file
# ds5day_file=$run_dir/dpscream_rce_large_3km_${run}_h0_last5days.nc
tracer_file=$run_dir/dpscream_rce_large_3km_${run}

# cdo -seltimestep,-20/-1 $ds_file $ds5day_file

cdo -fldmean -timmean -seltimestep,-20/-1 -selvar,Z3 $ds_file $run_dir/dpscream_rce_large_3km_${run}_Z3_mean.nc

tracer=BCU
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
cdo -mulc,-1 -ln -selvar,$tracer $ds5day_file $tracer_file # -seltimestep,-20/-1 

tracer=NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
cdo -mulc,-1 -ln -selvar,$tracer $ds5day_file $tracer_file # -seltimestep,-20/-1

tracer=W_NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
cdo -setname,$tracer -div -selvar,$tracer $ds5day_file -selvar,NUC $ds5day_file ${tracer_file}

tracer=NI_NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_hrs.nc
cdo -setname,$tracer -div -selvar,$tracer $ds5day_file -selvar,NUC $ds5day_file ${tracer_file}

echo "done"
