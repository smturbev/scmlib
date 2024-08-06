#!/bin/bash
#PBS -N zPreFSUT
#PBS -A UWAS0108
#PBS -l walltime=00:15:00
#PBS -q develop
#PBS -j oe
#PBS -k eod
#PBS -m be
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=8

module load cdo


# plot types can be the following:
# 'binned_by_crh', ...
plot_type="binned_by_crh"
ndays=5
# var = 3D: IWC NCU BCU W 
# var = 2D: LWCF SWCF TMQ FLUT FSUTOA LHFLX SHFLX  OMEGA500 

if [ $plot_type = 'binned_by_crh' ] ; then
    var="FSUTOA"
    sw_rad=.true.
    run="aa_default"
    declare -a BinArray=(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95)

    run_dir=/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/dpscream_rce_large_3km_$run/run
    crh_file=$run_dir/dpscream_rce_large_3km_${run}_TMQ_percs_2d_new.nc
    # last${ndays}days.nc
    if [ $var = 'W' ] ; then
        var_file=$run_dir/dpscream_rce_large_3km_${run}_W_last${ndays}days.nc
    else
        var_file=$run_dir/dpscream_rce_large_3km_${run}.eam.h0.2000-01-01-00000.nc
    fi 
    mask_temp=$run_dir/temp_${var}_mask_bins_crh.nc
    if [ $sw_rad ] ; then
        out_temp=$run_dir/temp_${var}_out_crh_bins_daytime
    else
        out_temp=$run_dir/temp_${var}_out_crh_bins
    fi

    for i in "${BinArray[@]}"; do
        echo $i" to "$(($i + 5))"for "$ndays" days"
        cdo -add -gec,$i $crh_file -ltc,$(($i + 5)) -selvar,"TMQ_percs" $crh_file $mask_temp
        # cdo -add -gec,$i $crh_file -ltc,$(($i + 5)) -mulc,100 -selvar,"TMQ_percs" $crh_file $mask_temp
        # cdo -timmean -fldmean -ifthen $mask_temp -seldate,2000-02-15T06:00:00,2000-02-19T18:00:00 -selvar,$var $var_file ${out_temp}_${i}.nc
        # cdo -timmean -fldmean -ifthen $mask_temp -seltimestep,$((-4*$ndays))/-2 -selvar,$var $var_file ${out_temp}_$(printf %02d $i).nc
        if [ $sw_rad ] ; then
            ## for SW radiation:
            cdo -timmean -fldmean -selhour,10,11,12,13,14 -ifthen $mask_temp -seltimestep,$((-4*$ndays))/-2 -selvar,$var $var_file ${out_temp}_$(printf %02d $i).nc
        else
            cdo -timmean -fldmean -ifthen $mask_temp -seltimestep,$((-4*$ndays))/-2 -selvar,$var $var_file ${out_temp}_$(printf %02d $i).nc
        fi
        
    done
fi


echo "done"
