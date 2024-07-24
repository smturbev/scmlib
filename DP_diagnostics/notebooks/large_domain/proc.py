''' proc.py
    created by Sami Turbeville
    on 7/22/2024

    records/holds processing tools used for scream_large runs
'''
# from matplotlib import cm
# import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
# from scipy import stats


def save_crh_percentiles(run_file, ndays=5):
    """ Returns the column humidity percentiles using TMQ
        (precipitable water or column-integrated water vapor).
    Inputs:
        - run_file  (str) is for the 3D file
    Outputs:
        - None, saves the TMQ_percentiles in a new file in 
            same directory as given file
    """
    var="IWC"
    chunks={'lev':64,'ncol':2000}
    bins = np.arange(0,1.01,0.05)
    bin_mids = (bins[1:]+bins[:-1])/2

    ds = xr.open_dataset(run_file, chunks=chunks).isel(time=slice(-ndays*4,-1))
    ds = ds[var]
    file2d = run_file.replace("h0","h1")
    crh = xr.open_dataset(file2d)["TMQ"].sel(time=ds.time)
    crh_perc_bins = crh.quantile(bins).values
    crh_percs = np.zeros(crh.shape)
    for i in range(len(bins)-1):
        print(i, bin_mids[i], crh_perc_bins[i], crh_perc_bins[i+1])
        crh_percs = np.where((crh>=crh_perc_bins[i]) & (crh<crh_perc_bins[i+1]),
                             bin_mids[i], crh_percs)
    crh_percs = np.repeat(crh_percs[:,np.newaxis,:], 128, axis=1)
    print("tmq_percs shape vs 3d shape", crh_percs.shape, ds.shape)
    crh_percs = xr.DataArray(crh_percs, dims=ds.dims, coords=ds.coords,
                             attrs={"name":"column-integrated water vapor percentiles"})
    crh_percs = xr.Dataset({"TMQ_percs":crh_percs}, attrs=crh_percs.attrs)
    run_dir, file_name = run_file.split(run_file.split("/")[-1])
    crh_percs.to_netcdf(run_dir+file_name.split(".eam")[0]+"TMQ_percs.nc")
    print("saved as "+run_dir+file_name.split(".eam")[0]+"TMQ_percs.nc")
    return

def cat_var_binned_by_crh(var, run_dir, file_base_name, new_file_name=None):
    """ Returns a concatenauted file of the binned variable 
        in a single file with bin edges included.

        Input:
            - run_dir (str): directory of files
            - file_base (str): base of the file (before the bin # and .nc)
            - new_file (str): new file name
        Output:
            - None: saves new file
    """
    print('getting', var, 'from', run_dir+file_base_name+"_i.nc")
    bins = np.arange(0,101,5)
    new_var = np.zeros((len(bins)-1,128))
    ds = xr.open_dataset(run_dir+file_base_name+"_0.nc")[var]
    print(ds.shape)
    for i in range(len(bins[:-1])):
        print(i, bins[i])
        new_var[i,:] = xr.open_dataset(run_dir+file_base_name+"_"+str(bins[i])+".nc")[var].isel(time=0).isel(x=0).isel(y=0).values
    new_var = xr.DataArray(new_var, dims=['bins','lev'], 
                           coords={'bins':((bins[1:]+bins[:-1])/2), 'lev':ds.lev})
    new_ds = xr.Dataset({var:new_var}).assign_coords({"bin_edges":bins})
    if new_file_name is not None:
        new_ds.to_netcdf(run_dir+new_file_name+".nc")
        print("saved as "+run_dir+new_file_name+".nc")
    return new_ds
        
    