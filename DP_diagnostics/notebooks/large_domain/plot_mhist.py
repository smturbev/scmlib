import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from matplotlib import colors
import util

# from dask.diagnostics import ProgressBar
# pbar = ProgressBar()
# pbar.register()

chunks = {"x":2000, "ncol":2000}
plt.rcParams.update({'font.size': 18})  # all to 18 if not specified other
qsmall=1e-8

####################################################################
### user input 
run = "dpscream_rce_large_3km_default304k_branch"
stat='norm'
### end user input
####################################################################

run_dir = f"/glade/derecho/scratch/sturbeville/DPSCREAM_simulations/{run}/run/"
file = run_dir + f"{run}_h0_last5days.nc"
print(file)

ds = xr.open_dataset(file, chunks=chunks, engine="netcdf4")[["T","CLDICE","NUMICE","Q"]]
print("getting t<-40degC and cldice>qsmall...")
ds = ds.where((ds["T"]<233.15)&(ds.CLDICE>qsmall))
print(ds.Q.shape)

print("calc Rice...")
x_array = util.calc_rice(ds.CLDICE, ds.NUMICE)
print(x_array.shape)

print("calc NI...")
y_array = util.calc_ni(ds.NUMICE, ds.Q, ds.lev*100, ds["T"])
print(y_array.shape)

del ds

print("flattening...")
x_array = x_array.values.flatten()
print("x done...")

n = len(x_array)
y_array = y_array.values.flatten()/1e6  # convert to cm-3
print("y done...")

xbins=np.linspace(0,100,100)
ybins=np.logspace(-5,2,70)

hist, xbins, ybins, _ = stats.binned_statistic_2d(x_array, y_array, None,
                                                  statistic='count', bins=[xbins,ybins])
if stat=="norm":
    print("normalized histogram...")
    hist = np.where(hist>0,hist,np.nan)
    hist = hist/n
print("hist sum", np.nansum(hist),stat)

xbins = (xbins[1:]+xbins[:-1])/2
ybins = (ybins[1:]+ybins[:-1])/2

print("plotting... ", end="")
fig, ax = plt.subplots(1, 1, figsize=(7,7), constrained_layout=True)
cf = ax.pcolormesh(xbins, ybins, (hist*100).T, cmap="magma_r",
                   shading='auto')
ax.set(yscale='log')
ax.set(xlim=[0,100], ylim=[5e-5,10])
ax.set(ylabel="ICNC (#/cm$^3$)", xlabel="Rice (um)")
plt.colorbar(cf, ax=ax, label=f"pdf", location="bottom", shrink=0.8, extend='max')
plt.savefig(f"../plots/large/micro_hist_run_{run}_{stat}.pdf")
plt.close()
