#!/home/disk/p/smturbev/miniconda3/envs/dyamond/bin/python

""" plot_tracers_cdf.py """


import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

plot_insitu_only=True  # plot it with ONLY the insitu ciruss cdf line
plot_wo_insitu=False  # plot it without the insitu cirrus line

if (plot_insitu_only & plot_wo_insitu):
    raise Exception('plot_insitu_only and plot_wo_insitu cannot both be True')

run = 'dpscream_rce_large_3km_d_default304k_wbranch'
file_dir = '/home/disk/olympus5/smturbev/DPSCREAM_simulations/ice_sensitivity_3d_variables/'
ice_file = file_dir + f'{run}_h0_last5days.nc'
bcu_file = f'{file_dir}{run}_BCU_hrs.nc'
nuc_file = f'{file_dir}{run}_NUC_hrs.nc'
# wnuc_file = f'{file_dir}{run}_W_NUC_hrs.nc'

color_dict = {'default':'gray', 'lsascent':'darkcyan', 'lpfrz':'cornflowerblue','lpls':'darkorange'}

nt = 10  # last nt timesteps
chunks = {'ncol':2000}
qsmall= 1e-8
plt.rc('font', size=16)

print(f'Loading files for {run}...')
print('... NUC...')
nuc_ds = xr.open_dataset(nuc_file, chunks=chunks).NUC.isel(time=slice(-nt,-1))
print('... BCU...')
bcu_ds = xr.open_dataset(bcu_file, chunks=chunks).BCU.isel(time=slice(-nt,-1))
print('... other 3d variables...')
ice_ds = xr.open_dataset(ice_file, chunks=chunks).IWC.isel(time=slice(-nt,-1))
# print('... W_NUC...')
# wnuc_ds = xr.open_dataset(wnuc_file, chunks=chunks).W_NUC.isel(time=slice(-nt,-1))

print('... assert the files have the same time dimension...')
assert nuc_ds.time[0] == bcu_ds.time[0]
assert bcu_ds.time[0] == ice_ds.time[0]
print('... done. Get cloudy points...')

nuc_ds = nuc_ds.where(ice_ds>qsmall)
bcu_ds = bcu_ds.where(ice_ds>qsmall)
n = np.nansum(ice_ds.values>qsmall)
print(n, 'total cloudy grid points with ice')
mask = ((bcu_ds>24)&(nuc_ds<bcu_ds-5))
ins_ci = np.where(mask, nuc_ds, np.nan)
n_ins = np.nansum(mask)
print(n_ins, 'total insitu cirrus grid points')

bins = np.arange(0,48)
bin_mids = (bins[1:]+bins[:-1])/2
print('Get histogram...')
if not(plot_insitu_only):
    print('... NUC hist...')
    nuc_hist, _ = np.histogram(nuc_ds.values.flatten(), bins=bins)
    print('... BCU hist...')
    bcu_hist, _ = np.histogram(bcu_ds.values.flatten(), bins=bins)
print('... Insitu cirrus hist...')
ins_hist, _ = np.histogram(ins_ci.flatten(), bins=bins)

print('Turn into cdf (cumsum)...')
if not(plot_insitu_only):
    print('... NUC cdf...')
    nuc_cdf = np.cumsum(nuc_hist)/n
    print('... BCU cdf...')
    bcu_cdf = np.cumsum(bcu_hist)/n
print('... Insitu cirrus cdf...')
ins_cdf = np.cumsum(ins_hist)/n_ins

print('print data points in cdf...')
if not(plot_insitu_only): print('nuc...',dict(zip(bin_mids, nuc_cdf)))
if not(plot_insitu_only): print('bcu...',dict(zip(bin_mids, bcu_cdf)))
print('...ins...',dict(zip(bin_mids, ins_cdf)))

print('Plotting...')
run_label=run.split('_')[-2]
print(run_label[:-4])
fig, ax = plt.subplots(1,1, figsize=(5,5))
if not(plot_insitu_only):
    ax.plot(bin_mids, nuc_cdf, label='time since nucleation')
    ax.plot(bin_mids, bcu_cdf, label='tiime since convection')
if not(plot_wo_insitu):
    if plot_insitu_only:
        if '300' in run_label:
            ls='solid'
        else:
            ls='dashed'
        ax.plot(bin_mids, ins_cdf, label=f'{run_label}',
                color=color_dict[run_label[:-4]], linestyle=ls)
        ax.set(title='insitu cirrus')
    else:
        ax.plot(bins_mids, ins_cdf, label='insitu ci. time since nuc')
ax.set(xlabel='hours', ylabel='cdf', ylim=[0,1], xlim=[0,48])
ax.grid(True)
ax.legend(loc=4)
if plot_insitu_only:
    savename=f'../plots/large/tracers_cdf_insitu_only_{run}.png'
else:
    savename=f'../plots/large/tracers_cdf_{run}.png'
plt.savefig(savename, dpi=150, transparent=True, bbox_inches='tight', pad_inches=0.2)
plt.close()
print(f'... Saved as {savename}\nDone.')