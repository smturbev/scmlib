#!/bin/tcsh
#PBS -N zDPSCREAM_SMALL_test_bcu_aa
#PBS -A UWAS0108
#PBS -l walltime=00:40:00
#PBS -q develop
#PBS -j oe
#PBS -k eod
#PBS -m e
#PBS -M smturbev@uw.edu
#PBS -l select=1:ncpus=8:mpiprocs=8

# This script compiles scream and sets a DP-SCREAM case running
#   based on the RCE setup below.

# Set up modules
source /glade/u/apps/derecho/23.09/spack/opt/spack/lmod/8.7.24/gcc/7.5.0/c645/lmod/lmod/init/tcsh
module purge
module load ncarenv/23.09 craype intel/2024.0.2 ncarcompilers cray-mpich hdf5 netcdf cmake parallel-netcdf parallelio/2.6.2-static
#######################################################################
#######################################################################
#######  Script to run SCREAM in doubly periodic (DP) mode
#######  RCE
#######  Radiative Convective Equilibrium
#######
#######  Script Author: P. Bogenschutz (bogenschutz1@llnl.gov)

# Modified for derecho-cpu by P. Blossey (pblossey@uw.edu)
#######################################################
#######  BEGIN USER DEFINED SETTINGS

  # Set the name of your case here
  setenv casename scream_dp_RCE_SMALL_v0 #dpscream_rce_small_3km_aa_default
  
# Set the case directory here
  setenv casedirectory /glade/derecho/scratch/$USER/DPSCREAM_simulations

  # Directory where inputdata lives
  setenv inputdata_dir /glade/work/$USER/E3SM/inputdata

  # Directory where code lives
  setenv code_dir /glade/u/home/$USER

  # Code tag name
  setenv code_tag scream

  # Name of machine you are running on (i.e. cori, anvil, etc)
  setenv machine derecho-cpu

  # Name of project to run on, if submitting to queue
  setenv projectname UWAS0108

  # Name of queue for job submission
  setenv job_queue main
  setenv job_priority economy

  # set email for error messages
  set email = smturbev@uw.edu

  # Set to debug queue?
  # - Some cases are small enough to run on debug queues
  # - Setting to true only supported for NERSC and Livermore Computing,
  #   else user will need to modify script to submit to debug queue
  setenv debug_queue false

  # Set number of processors to use, should be less than or equal
  #   to the total number of elements in your domain.
  set num_procs = 128 # 512

  # set walltime
  set walltime = '12:00:00'

  ## SET DOMAIN SIZE AND RESOLUTION:
  # - Note that these scripts are set to run with dx=dy=3.33 km
  # which is the default SCREAM resolution.

  # To estimate dx (analogous for dy):
  # dx = domain_size_x / (num_ne_x * 3)
  # (there are 3x3 unique columns per element, hence the "3" factor)

  # Set number of elements in the x&y directions
  set num_ne_x = 12
  set num_ne_y = 12

  # Set domain length [m] in x&y direction
  set domain_size_x = 120000
  set domain_size_y = 120000

  # BELOW SETS RESOLUTION DEPENDENT SETTINGS
  # (Note that all default values below are appropriate for dx=dy=3.33 km and do not
  #  need to be modified if you are not changing the resolution)

  # SET MODEL TIME STEPS
  #  -NOTE that if you change the model resolution,
  #  it is likely the model and physics time steps will need to be adjusted.
  #  As a rule, a factor of 2 increase in resolution should equate to a factor of 2
  #  decrease of the model time steps.

  # model and physics time step [s]
  set model_dtime = 100

  # dynamics time step [s]
  #  should divide evenly into model_dtime
  set dyn_dtime = 8.333333333333333d0

  # SET SECOND ORDER VISCOSITY NEAR MODEL TOP
  #  NOTE that if you decrease resolution you will also need to reduce
  #  the value of "nu_top" (second-order viscosity applied only near model top).
  #  Rule of thumb is that a factor of 2 increase in resolution should equate to a
  #  factor of 2 decrease for this value

  # second order visocosity near model top [m2/s]
  set nu_top_dyn = 1.e4

####### END (mandatory) USER DEFINED SETTINGS, but...
####### Likely POSSIBLE EXCEPTIONS (not limited to):
#######  - If the user wants to add additional output, for example, the EAM
#######	   namelist (user_nl_eam) should be modified below to accomodate for this.
#######  - The user wants to run a subset of the selected case.
#######
#######  - NOTE ON DEFAULT OUTPUT
#######    - *eam.h0* tapes contain the the default output averaged daily
#######      (for multi-day cases) or hourly (for shorter boundary layer
#######      cloud cases)
#######    - *eam.h1* tapes contain instantaneous 2D fields output hourly
#######    - ALL/any of this can be modified by the user based on needs
###########################################################################
###########################################################################
###########################################################################

# Case specific information kept here
  set lat = 0.0 # latitude
  set lon = 0.0 # longitude
  set do_iop_srf_prop = .false. # Use surface fluxes in IOP file?
  set do_iop_nudge_tq = .false. # Relax T&Q to observations?
  set do_iop_nudge_uv = .false. # Relax U&V to observations?
  set do_iop_subsidence = .false. # compute LS vertical transport?
  set do_turnoff_swrad = .false. # Turn off SW calculation
  set do_turnoff_lwrad = .false. # Turn off LW calculation
  set startdate = 2000-01-01 # Start date in IOP file
  set start_in_sec = 0 # start time in seconds in IOP file
  set stop_option = ndays
  set stop_n = 5
  set iop_file = RCE_iopfile_4scam_no-mean_ascent.nc #IOP file name
  set sst_val = 300 # set constant SST value (ONLY valid for RCE case)
  set p3_new_icenuc = .false. # Turn off new ice nucleation scheme in P3
# End Case specific stuff here

  # Location of IOP file
  set iop_path = atm/cam/scam/iop

  set PROJECT=$projectname
  set E3SMROOT=${code_dir}/${code_tag}

  cd $E3SMROOT/cime/scripts

  set compset=F2010-SCREAM-HR

  # Note that in DP-SCREAM the grid is set ONLY to initialize
  #  the model from these files
  set grid=ne30_ne30

  set CASEID=$casename

  set CASEDIR=${casedirectory}/$CASEID

  set run_root_dir = $CASEDIR
  set temp_case_scripts_dir = $run_root_dir/case_scripts

  set case_scripts_dir = $run_root_dir/case_scripts
  set case_build_dir   = $run_root_dir/build
  set case_run_dir     = $run_root_dir/run

# Create new case
  ./create_newcase -case $casename --script-root $temp_case_scripts_dir -mach $machine -project $PROJECT -compset $compset -res $grid
  cd $temp_case_scripts_dir

  ./xmlchange JOB_WALLCLOCK_TIME=$walltime

# Define executable and run directories
  ./xmlchange --id EXEROOT --val "${case_build_dir}"
  ./xmlchange --id RUNDIR --val "${case_run_dir}"
  ./xmlchange --id DIN_LOC_ROOT --val "${inputdata_dir}"
  ./xmlchange --id JOB_QUEUE --val "${job_queue}"
  ./xmlchange --id JOB_PRIORITY --val "${job_priority}"

# Set to debug, only on certain machines
  if ($debug_queue == 'true') then
    if ($machine =~ 'cori*') then
      ./xmlchange --id JOB_QUEUE --val 'debug'
    endif

    if ($machine == 'quartz' || $machine == 'syrah') then
      ./xmlchange --id JOB_QUEUE --val 'pdebug'
    endif
  endif

# Get local input data directory path
  set input_data_dir = `./xmlquery DIN_LOC_ROOT --value`

# need to use single thread
  set npes = $num_procs
  foreach component ( ATM LND ICE OCN CPL GLC ROF WAV )
    ./xmlchange  NTASKS_$component=$npes,NTHRDS_$component=1,ROOTPE_$component=0
  end

# CAM configure options.  Set to SCREAM default settings.
  set CAM_CONFIG_OPTS="-phys default -scam -dpcrm_mode -nlev 128 -shoc_sgs -microphys p3 -rad rrtmgp -chem spa -cldera_passive_trcs"

  ./xmlchange CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS"
  ./xmlchange --id CAM_CONFIG_OPTS --append --val='-rce -aquaplanet'

  # Always run with the theta-l version of HOMME, the default for SCREAM
  ./xmlchange CAM_TARGET=theta-l

# if we want to turn off SW radiation, then set appropriate namelist settings here
  if ($do_turnoff_swrad == '.true.') then
    set iradsw_in = 0
  else
    set iradsw_in = 3
  endif

# if we want to turn off LW radiation, then set appropriate namelist settings here
  if ($do_turnoff_lwrad == '.true.') then
    set iradlw_in = 0
  else
    set iradlw_in = 3
  endif

# Compute maximum allowable number for processes (number of elements)
  set dyn_pes_nxny = `expr $num_ne_x \* $num_ne_y`

# Runtime specific namelist information
cat <<EOF >> user_nl_eam
 use_gw_front = .false.
 use_gw_oro = .false.
 use_gw_convect = .false.
 deep_scheme = 'off'
 convproc_do_aer = .false.
 iop_dosubsidence = $do_iop_subsidence
 iop_nudge_tq = $do_iop_nudge_tq
 iop_nudge_uv = $do_iop_nudge_uv
 history_aerosol = .false.
 micro_tend_output = .true.
 theta_hydrostatic_mode = .false.
 tstep_type = 9
 do_new_bg_lp_frz = $p3_new_icenuc
 dep_scaling_small = 1.0
 sed_scaling_small = 1.0
 scale_all_ice = .false.
fexcl1='BRUNT','FICE','EXTINCT','FREQI','FREQL','FREQR','FREQS','RELVAR','UU','VQ','VT','VU','VV','AODABS','AODABSBC','AODALL','AODBC','AODDUST','AODDUST1','AODDUST3','AODMODE1','AODMODE2','AODMODE3','AODNIR','AODPOM','AODSO4','AODSOA','AODSS','AODUV','AODVIS','BURDEN1','BURDEN2','BURDEN3','CCN3' fincl2='CAPE','CIN','CLDLOW','CLDMED','CLDHGH','CLDTOT','CDNUMC','DTENDTH','DTENDTQ','FLDS','FLNS','FLNSC','FLNT','FLNTC','FLUT','FLUTC','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTC','FSNTOA','FSNTOAC','FSUTOA','FSUTOAC','LHFLX','SHFLX','LWCF','SWCF','OMEGA500','PRECC','PRECL','PS','QREFHT','SOLIN','TAUX','TAUY','TGCLDCWP','TGCLDIWP','TGCLDLWP','TH7001000','TMQ','TREFHT','TS','WINDSPD_10M','crm_grid_x','crm_grid_y'
 fincl1='OMEGA','DYN_OMEGA','QRL','QRS','CLDICE','WSUB'
 mfilt = 5000, 5000
 nhtfrq = -6, -6
 avgflag_pertape='I','I'
 scmlat = $lat
 scmlon = $lon
 iradsw = $iradsw_in
 iradlw = $iradlw_in
 scm_iop_srf_prop = $do_iop_srf_prop
 iopfile = '$input_data_dir/$iop_path/$iop_file'
 pertlim = 0.001
 iop_perturb_high = 900.0D0
 ncdata='$input_data_dir/atm/cam/inic/homme/cami_mam3_Linoz_ne30np4_SCREAM_L128_c160214.nc'
EOF

# Timestepping stuff related to DP-SCREAM
# NOTE, if you change resolution from default it may be required to
#  change some of these settings.
cat <<EOF >> user_nl_eam
 dyn_npes=$dyn_pes_nxny
 se_tstep=$dyn_dtime
 cldfrc_iceopt = 7
 transport_alg = 0
 dt_tracer_factor = 1
 hypervis_subcycle_q = 1
 dt_remap_factor = 1
 nu             =   0.216784
 nu_top         =  $nu_top_dyn
 se_ftype       = 2
 cubed_sphere_map = 2
 cld_macmic_num_steps =  1
 shoc_timestep = -1
EOF

# Settings related to domain size and resolution
cat <<EOF >> user_nl_eam
 mesh_file = 'none'
 se_ne_x = $num_ne_x
 se_ne_y = $num_ne_y
 se_lx = $domain_size_x
 se_ly = $domain_size_y
EOF

# avoid the monthly cice file from writing as this
#   appears to be currently broken for SCM
cat <<EOF >> user_nl_cice
  histfreq='y','x','x','x','x'
EOF

# Turn on UofA surface flux scheme
cat <<EOF>> user_nl_cpl
  ocn_surface_flux_scheme = 2
EOF


# Modify the run start and duration parameters for the desired case
 ./xmlchange RUN_STARTDATE="$startdate",START_TOD="$start_in_sec",STOP_OPTION="$stop_option",STOP_N="$stop_n"
# For a branched run, input the reference directory here
 # ./xmlchange RUN_TYPE="branch",RUN_REFCASE=$branch_case,RUN_REFDATE="2000-01-06",RUN_REFDIR="/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/${branch_case}/run/"

# Compute number of columns needed for component model initialization
  set comp_mods_nx = `expr $num_ne_x \* $num_ne_y \* 9`

# Modify the latitude and longitude for the particular case
  ./xmlchange PTS_MULTCOLS_MODE="TRUE",PTS_MODE="TRUE",PTS_LAT="$lat",PTS_LON="$lon"
  ./xmlchange MASK_GRID="USGS",PTS_NX="${comp_mods_nx}",PTS_NY=1
  ./xmlchange ICE_NX="${comp_mods_nx}",ICE_NY=1

  ./xmlchange DOCN_AQPCONST_VALUE=$sst_val
  ./xmlchange DOCN_MODE="sst_aquap_constant"


# Set model timesteps

  @ ncpl = 86400 / $model_dtime
  ./xmlchange ATM_NCPL=$ncpl
  ./xmlchange CAM_NAMELIST_OPTS="dtime=$model_dtime"

  ./case.setup

# Write restart files at the end of model simulation
  ./xmlchange PIO_TYPENAME="pnetcdf"
  ./xmlchange REST_OPTION="end"

# Build the case
  ./case.build

# Submit the case
  ./case.submit --mail-user $email -M begin,end

  exit
