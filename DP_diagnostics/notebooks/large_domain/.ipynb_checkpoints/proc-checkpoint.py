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


def save_crh_percentiles(run_file, ndays=5, save2d=False):
    """ Returns the column humidity percentiles using TMQ
        (precipitable water or column-integrated water vapor).
    Inputs:
        - run_file  (str) is for the 3D file
    Outputs:
        - None, saves the TMQ_percentiles in a new file in 
            same directory as given file
    """
    chunks={'lev':64,'ncol':2000}
    bin_edges = np.arange(0,101,5)
    bin_mids = (bin_edges[1:]+bin_edges[:-1])/2

    ds = xr.open_dataset(run_file, chunks=chunks).isel(time=slice(-ndays*4,-1))[["TMQ"]]
    ds_time = ds.TMQ.time
    ds_dims = ds.TMQ.dims
    ds_coords = ds.TMQ.coords

    crh = ds["TMQ"].chunk(dict(ncol=-1))

    del ds

    crh_perc_bins = crh.quantile(bin_edges/100).values[1:-1]
    crh_percs = np.zeros(crh.shape)

    crh_percs = np.where((crh<crh_perc_bins[0]),
                         bin_mids[0], crh_percs)
    print("crh < 5% i.e.,", bin_edges[0],"to", bin_edges[1],
          "with a bin mid of", bin_mids[0], "less than", crh_perc_bins[0])
    print("... ",(crh<crh_perc_bins[0]).sum().values)

    for i in range(1,len(bin_edges)-2):
        print(i, "crh is between", int(bin_edges[i]),"and", int(bin_edges[i+1]),
              "with percentiles between", int(crh_perc_bins[i-1]), int(crh_perc_bins[i]),
              "with a bin mid of", (bin_mids[i]))
        crh_percs = np.where((crh>=crh_perc_bins[i-1]) & (crh<crh_perc_bins[i]),
                             bin_mids[i], crh_percs)
        print("... ", ((crh>=crh_perc_bins[i-1]) & (crh<crh_perc_bins[i])).sum().values)

    crh_percs = np.where((crh>=crh_perc_bins[-1]),
                             bin_mids[-1], crh_percs)
    print("crh >= 95% i.e.,", bin_edges[-2],"to", bin_edges[-1],
          "with a bin mid of", bin_mids[-1])
    print("... ", (crh>=crh_perc_bins[-1]).sum().values)

    if save2d:
        crh_percs = xr.DataArray(crh_percs, dims=ds_dims, coords=ds_coords,
                         attrs={"name":"column-integrated water vapor percentiles"})
        savename="_2d"
    else:
        crh_percs = np.repeat(crh_percs[:,np.newaxis,:], 128, axis=1)
        print("tmq_percs shape vs 3d shape", crh_percs.shape, ds_dims, ds_coords)
        crh_percs = xr.DataArray(crh_percs, dims=ds_dims, coords=ds_coords,
                                 attrs={"name":"column-integrated water vapor percentiles"})
        savename=""
    crh_percs = xr.Dataset({"TMQ_percs":crh_percs}, attrs=crh_percs.attrs)

    base_name = run_file.split("/")[-1]
    run_dir = run_file.split(base_name)[0]
    save_name = run_dir+base_name.split(".eam")[0]+f"_TMQ_percs{savename}_newnewnew.nc"
    print("saved as "+save_name)
    crh_percs.to_netcdf(save_name)
    return crh_percs

def cat_var_binned_by_crh(var, run_dir, file_base_name, new_file_name=None, save2d=False):
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
    # ds = xr.open_dataset(run_dir+file_base_name+"_0.nc")[var]
    # print(ds.shape)
    for i in range(len(bins[:-1])):
        file_name = run_dir+file_base_name+f"_{bins[i]:02d}.nc"
        new_var[i,:] = xr.open_dataset(file_name)[var][0,0,0].values
        print(i, bins[i], file_name, new_var[i,0])
    if save2d:
        new_var = xr.DataArray(new_var[:,0], dims=['bins'], 
                               coords={'bins':((bins[1:]+bins[:-1])/2)})
    else:
        new_var = xr.DataArray(new_var, dims=['bins','lev'], 
                               coords={'bins':((bins[1:]+bins[:-1])/2), 'lev':new_var.lev})
    new_ds = xr.Dataset({var:new_var}).assign_coords({"bin_edges":bins})
    if new_file_name is not None:
        new_ds.to_netcdf(run_dir+new_file_name+".nc")
        print("saved as "+run_dir+new_file_name+".nc")
    return new_ds
        
    