#!/bin/bash -v
#PBS -N ztracLSascent300K
#PBS -A UWAS0108
#PBS -l walltime=01:00:00
#PBS -q develop
#PBS -j oe
#PBS -k eod
#PBS -m be
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=87

module load cdo
 
run=CoopLS_300K
run_dir=/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/dpscream_rce_large_3km_${run}/run/
ds_file=$run_dir/dpscream_rce_large_3km_${run}.eam.h0.2000-01-01-00000.nc
ds5day_file=$run_dir/dpscream_rce_large_3km_${run}_h0_days40-50.nc

# cdo -seltimestep,-20/-1 $ds_file $ds5day_file
cdo -seldate,2000-02-09,2000-02-19 -selvar,NUC,BCU,W_NUC,NI_NUC,IWC,CLDICE,CLDLIQ,NUMICE,T,OMEGA,Q,Z3,QRS,QRL,TMQ,TGCLDIWP,TGCLDLWP,LWCF,SWCF,crm_grid_x,crm_grid_y $ds_file $ds5day_file

cdo -fldmean -timmean -selvar,Z3 $ds5day_file $run_dir/dpscream_rce_large_3km_${run}_Z3_mean.nc

tracer=BCU
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_days40-50.nc
cdo -mulc,-1 -ln -selvar,$tracer $ds5day_file $tracer_file # -seltimestep,-20/-1 

tracer=NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_days40-50.nc
cdo -mulc,-1 -ln -selvar,$tracer $ds5day_file $tracer_file # -seltimestep,-20/-1

tracer=W_NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_days40-50.nc
cdo -setname,$tracer -div -selvar,$tracer $ds5day_file -selvar,NUC $ds5day_file ${tracer_file}

tracer=NI_NUC
tracer_file=$run_dir/dpscream_rce_large_3km_${run}_${tracer}_days40-50.nc
cdo -setname,$tracer -div -selvar,$tracer $ds5day_file -selvar,NUC $ds5day_file ${tracer_file}

echo "done"
