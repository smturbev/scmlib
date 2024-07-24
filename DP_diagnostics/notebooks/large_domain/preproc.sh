#!/bin/bash
#PBS -N zpreproc_cdo
#PBS -A UWAS0108
#PBS -l walltime=00:15:00
#PBS -q main
#PBS -j oe
#PBS -k eod
#PBS -m be
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=8

module purge
module load cdo

# plot types can be the following:
# 'binned_by_crh', ...
plot_type="binned_by_crh"

if [ $plot_type = 'binned_by_crh' ] ; then
    var="W_NUC"
    declare -a BinArray=(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95)

    run_dir=/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/dpscream_rce_large_3km_aa_default/run
    crh_file=$run_dir/dpscream_rce_large_3km_aa_default_TMQ_percs.nc
    var_file=$run_dir/dpscream_rce_large_default.eam.h0.2000-01-01-00000.nc
    mask_temp=$run_dir/temp_mask_bins_crh.nc
    out_temp=$run_dir/temp_var_out_crh_bins

    for i in "${BinArray[@]}"; do
        echo $i" to "$(($i + 5))
        cdo -add -gec,$i $crh_file -ltc,$(($i + 5)) -mulc,100 -selvar,"TMQ_percs" $crh_file $mask_temp
        # if [ $var = 'NUC' ] ; then
        #     cdo -timmean -fldmean -ifthen $mask_temp -seldate,2000-02-15T06:00:00,2000-02-19T18:00:00 -selvar,$var $var_file ${out_temp}_${i}.nc
        # elif [ $var = 'BCU' ] ; then
        #     cdo -timmean -fldmean -ifthen $mask_temp -seldate,2000-02-15T06:00:00,2000-02-19T18:00:00 -selvar,$var $var_file ${out_temp}_${i}.nc
        # else
        cdo -timmean -fldmean -ifthen $mask_temp -seldate,2000-02-15T06:00:00,2000-02-19T18:00:00 -selvar,$var $var_file ${out_temp}_${i}.nc
        # fi
    done
fi

echo "done"
