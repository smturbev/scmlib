''' util_small.py
    created by Sami Turbeville
    on 9/21/23
'''
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy import stats

# universal constants
R = 287.058
G = 9.80665


def get_comp_names(comp_name):
    ''' Takes a string and returns the run codes, names, and colors'''
    if comp_name=="small":
        runs = ["f_default","g_halfsed","g_2xsed","g_halfdep", "g_2xdep"]
        run_names = ["Default","1/2x sed","2x sed","1/2x dep", "2x dep"]
        colors = ["darkred","lightskyblue","mediumblue","lightgreen","darkgreen"]
    elif comp_name=="all":
        runs = ["f_default","h_halfsed_all","h_2xsed_all","h_halfdep_all", "h_2xdep_all"]
        run_names = ["Default","1/2x sed (all)","2x sed (all)","1/2x dep (all)", "2x dep (all)"]
        colors = ["darkred","lightskyblue","mediumblue","lightgreen","darkgreen"]
    elif comp_name=="sed":
        runs = ["f_default","g_halfsed","g_2xsed","h_halfsed_all", "h_2xsed_all"]
        run_names = ["Default","1/2x sed (small)","2x sed (small)","1/2x sed (all)", "2x sed (all)"]
        colors = ["darkred","lightskyblue","mediumblue","lightgreen","darkgreen"]
    elif comp_name=="dep":
        runs = ["f_default","g_halfdep","g_2xdep","h_halfdep_all", "h_2xdep_all", "i_lsascent_a"]
        run_names = ["Default","1/2x dep (small)","2x dep (small)","1/2x dep (all)", "2x dep (all)", "LS ascent"]
        colors = ["darkred","lightskyblue","mediumblue","lightgreen","darkgreen", "darkviolet"]
    elif comp_name=="lsascent":
        runs = ["f_default", "i_lsascent_a", "i_halfomega"] # "c_lp2005"
        run_names = ["Default", "LS ascent", "Half ascent"]
        colors = ["darkred","darkviolet", "mediumvioletred"]
    elif comp_name=="nuc":
        runs = ["lpfrz_ac", "cimoh_lp_nohet"]  # , "i_bgicenuc"] # "c_lp2005"
        run_names = ["LP2005", "Ci Mohler"]  # , "LP + Ci Moh", "LP+ci"]
        colors = ["darkred","lightcoral", "magenta"]  # ,"darkorange","magenta"]
    elif comp_name=="sst":
        runs = ["f_default","j_304K","j_296K"]
        run_names = ["Default (300K)","Warm SST (304K)","Cool SST (296K)"]
        colors = ["grey","maroon","lightblue"]
    elif comp_name=="tracers_nuc":
        runs = ["cooper.nuc_tend_cooper_p3","cooper.nuc_tend_hom_p3",
                "lp05.nuc_tend_lp05_p3", "lp05.nuc_tend_hom_p3"]
        run_names= ["Cooper old", "Cooper new",
                    "LP05 old", "LP05 new"]
        colors=["lightblue", "darkblue", "lightgreen", "darkgreen"]
    else:
        raise Exception(f"comp_name {comp_name} not defined")
    return (runs, run_names, colors)


def omega2w(omega, p, t):
    ''' converts vertical velocity from pressure coords
        to m/s.

        input:
         - omega : array in units of Pa/s
         - p     : pressure array in Pa
         - t     : temperature array in K
        output:
         - w     : vertical velocity array in units of m/s
    '''
    rho = p / (R*t)
    w = -omega/(rho*G)
    return w


def w2omega(w, p, t):
    ''' converts vertical velocity from m/s
        to pressure coords (Pa/s).

        input:
         - w     : array in units of m/s
         - p     : pressure array in Pa
         - t     : temperature array in K
        output:
         - omega : vertical velocity array in units of Pa/s
    '''
    rho = p / (R*t)
    omega = -w*rho*G
    return omega


def calc_rice(qi, ni):
    rho = 920  # kg/m3
    # qi = qi  # kg/kg
    # ni = ni  # 1/kg
    # r_ice = np.where((ni>1e-5),(3*qi/(4*np.pi*rho*ni))**(1/3),0)*1e6  # um
    r_ice = (3*qi/(4*np.pi*rho*ni.where(ni>1e-5)))**(1/3)
    r_ice = r_ice * 1e6  # um
    return r_ice

def calc_ni(numice, qv, p, t):
    ni = numice * calc_rho(qv,p,t)  # 1/kg * kg/m3 = 1/m3
    return ni


def calc_rho(qv, p, t):
    """calculates density of air in kg/m3 for given input"""
    Tv = (1 + 0.61*qv) * t  # K
    rho = p / (R*Tv)  # kg/m3
    return rho


def calc_rhice(ds, varQ="Q", varT="T", z_units="hPa", magnus=True):
    """ input: xarray with variables Q and T
        output: xarray of rh wrt ice
    """
    if not(magnus):
        # from wv_sat_method.f90
        # (good down to 110 K)
        t = ds[varT]
        esi = np.exp(9.550426 - (5723.265 / t) + (3.53068 * np.log(t)) - (0.00728332 * t))
    else:
        t = ds[varT]-273.15
        esi = 611.21 * np.exp((22.587*t)/(t+273.86))
    if z_units=="hPa":
        z = ds.lev*100
    else:
        z=ds.lev
    wsi = (0.622 * esi) / z
    wi = ds[varQ] / (1 - ds[varQ])
    rh_ice = wi/wsi * 100
    return rh_ice

def calc_rhliq(ds, varQ="Q", varT="T", z_units="hPa"):
    """ input: xarray with variables Q and T
        output: xarray of rh wrt liquid
    """
    # # Magnus formula for saturation vapor pressure
    # t = ds[varT] - 273.15  # convert K to degC
    # e_s = 610.94 * np.exp((17.625*t)/(t+243.04))
    # ----
    # from EAM wv_sat_method.F90
    # (good for 123 < T < 332 K)
    t = ds[varT]
    es = ( np.exp(54.842763 - (6763.22 / t) - (4.210 * np.log(t)) +
           (0.000367 * t) + (np.tanh(0.0415 * (t - 218.8)) *
           (53.878 - (1331.22 / t) - (9.44523 * np.log(t)) +
           0.014025 * t)))
         )
    if z_units=="hPa":
        z = ds.lev*100
    else:
        z=ds.lev
    ws = (0.622 * es) / z
    w = ds[varQ] / (1 - ds[varQ])
    return( w/ws * 100 )

def calc_rhice_from_rhliq(ds, varT='T', varRH='RELHUM'):
    """ convert relative humidity wrt liquid to ice 
    using Magnus vapor pressure approximations"""
    t = ds[varT] # K
    es = (np.exp(54.842763 - (6763.22 / t) - (4.210 * np.log(t)) +
          (0.000367 * t) + (np.tanh(0.0415 * (t - 218.8)) *
          (53.878 - (1331.22 / t) - (9.44523 * np.log(t)) +
          0.014025 * t)))
          )
    esi = np.exp(9.550426 - (5723.265 / t) + (3.53068 * np.log(t)) - (0.00728332 * t))
    return (esi/es)*ds[varRH]

def plot_microhist(x_array, y_array, savename=None):
    """ plots the x vs y joint histogram (normalized pdf)

        Input: 
            x_array - Rice in um
            y_array - ICNC in 1/kg
            savemane - optional
        Output: 
            ax - figure axis
            saves figure as savename if specified
    """
    print(x_array.shape, y_array.shape)
    print("flattening...")
    x_array = x_array.values.flatten()
    y_array = y_array.values.flatten()/1e6  # convert to cm-3
    print(x_array.shape, y_array.shape)
    n = len(x_array)
    x_array = x_array[~np.isnan(y_array)]
    y_array = y_array[~np.isnan(y_array)]
    print(x_array.shape, y_array.shape)
    print(n)
    xbins=np.linspace(0,100,100)
    ybins=np.logspace(-5,2,70)
    
    hist, xbins, ybins, _ = stats.binned_statistic_2d(x_array, y_array, None,
                                                  statistic="count", bins=[xbins,ybins])
    hist = np.where(hist>0,hist,np.nan)
    hist = hist/n
    print(np.nansum(hist),"counts...",end="")
    xbins = (xbins[1:]+xbins[:-1])/2
    ybins = (ybins[1:]+ybins[:-1])/2
    
    print("plotting... ", end="")
    fig, ax = plt.subplots(1, 1, figsize=(6,6), constrained_layout=True)
    cf = ax.pcolormesh(xbins, ybins, (hist*100).T, cmap="magma_r", 
                      shading='auto',
                      # vmin = 0, vmax=0.2  #0.02, 0.3 std
                     )
    ax.set(yscale="log", xlabel="R$_{ice} (\mu$m)", 
       title=f"{savename} T<-40$^\circ$C",
       xlim=[0,100],ylim=[5e-5,10])
    ax.set_ylabel("ICNC (#/cm$^3$)")
    plt.colorbar(cf, ax=ax, label="%", location="bottom", shrink=0.8)
    if savename is not None:
        plt.savefig(f"../plots/small/micro_hist_run_{savename}.pdf")
    plt.show()
    print("done")
    return ax


def dennisplot(stat, olr, alb, var=None, xbins=None, ybins=None,
               levels=None, model="model", region="TWP", var_name=None, units="units",
               cmap=cm.ocean_r, ax=None, save=False, colorbar_on=True, fs=20):
    ''' Returns axis with contourf of olr and albedo.

    Parameters:
        - stat (str)   : - 'difference' returns contourf of the difference
                            between the first minus the second in the tuple
                         - 'density' returns density plot of olr-alb joint
                            histogram (pdf), or
                         - statistic for scipy.stats.binned_statistic_2d
        - olr (array)  : 1-D array of OLR values (from 85-310 W/m2),
        - alb (array)  : 1-D array of Albedo values (from 0-1),
        - var (array)  : 1-D array (optional if stat=density or difference)
        - colorbar_on (bool)
                       : returns a tuple of ax, mappable_countour if False

    Returns:
        - ax (plt.axis): axis with plot
        - cs (mappable): returned value from plt.contourf, if colorbar_on = False

    Note: Values for mean sw downward flux at toa from 
              http://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture11%20--%20Insolation.html. 
    '''
    if xbins is None:
        xbins = np.linspace(70,320,26)
    if ybins is None:
        ybins = np.linspace(0,0.8,33)
    if levels is None:
        if stat=="difference":
            levels = np.arange(-0.9,1,0.2)
        else:
            levels = np.arange(-3,-1.2,0.1)
    if stat=="difference":
        print("difference")
        olr0, olr1 = olr
        alb0, alb1 = alb
        olr0 = olr0[~np.isnan(alb0)]
        alb0 = alb0[~np.isnan(alb0)]
        alb0 = alb0[~np.isnan(olr0)]
        olr0 = olr0[~np.isnan(olr0)]
        olr1 = olr1[~np.isnan(alb1)]
        alb1 = alb1[~np.isnan(alb1)]
        alb1 = alb1[~np.isnan(olr1)]
        olr1 = olr1[~np.isnan(olr1)]
        hist0, xedges, yedges = np.histogram2d(olr0,alb0,bins=(xbins,ybins))
        nan_len = np.sum(~np.isnan(alb0))
        hist0 = hist0/nan_len
        print(nan_len)
        hist1, xedges, yedges = np.histogram2d(olr1,alb1,bins=(xbins,ybins))
        nan_len = np.sum(~np.isnan(alb1))
        hist1 = hist1/nan_len
        print(nan_len)
        binned_stat = hist0-hist1
    else:
        if stat=='density':
            if (olr.shape!=alb.shape):
                raise Exception("shapes don't match: olr %s, alb %s."%(olr.shape, alb.shape))
            olr = olr[~np.isnan(alb)]
            alb = alb[~np.isnan(alb)]
            alb = alb[~np.isnan(olr)]
            olr = olr[~np.isnan(olr)]
            # check for nans
            binned_stat, xedges, yedges = np.histogram2d(olr, alb, bins=(xbins,ybins))
            nan_len = xr.DataArray(alb).count().values
            binned_stat = binned_stat/(nan_len)
            print(nan_len)
        else:
            if (olr.shape!=var.shape):
                raise Exception("shapes don't match: olr %s, alb %s, %s %s."%(olr.shape, alb.shape, var_name, var.shape))
            var = var[~np.isnan(alb)]
            olr = olr[~np.isnan(alb)]
            alb = alb[~np.isnan(alb)]
            var = var[~np.isnan(olr)]
            alb = alb[~np.isnan(olr)]
            olr = olr[~np.isnan(olr)]
            alb = alb[~np.isnan(var)]
            olr = olr[~np.isnan(var)]
            var = var[~np.isnan(var)]
            var = var[~np.isnan(olr)]
            binned_stat, xedges, yedges, nbins = stats.binned_statistic_2d(olr, alb, var, bins=(xbins,ybins), statistic=stat)
    xbins2, ybins2 = (xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2
    if ax is None:
        ax = plt.gca()
    if stat=="difference":
        csn = ax.contourf(xbins2, ybins2, binned_stat.T*100, levels, cmap=cmap, extend='both')
    else:
        csn = ax.contourf(xbins2, ybins2, np.log10(binned_stat.T), levels, cmap=cmap, extend='both')
        ax.contour(csn, colors='k', linestyles='solid', linewidths=1)
    if region=="NAU":
        ax.plot([80,317],[0.57,0.],label="Neutral CRE", color='black')
        # calculated in turbeville et al., 2022
    elif region=="TWP":
        ax.plot([80,309],[0.55,0.],label="Neutral CRE", color='black')
        # calculated in turbeville et al., 2022
    else:
        ax.plot([80,320], [0.75,0.2], label="Neutral CRE", color='black')
        # calculated in turbeville et al., 2022
    ax.grid()
    ax.set_xticks([100,150,200,250,300])
    ax.set_ylim([0.05,0.8])
    ax.set_xlim([80,310])
    ax.set_xlabel('OLR(W m$^{-2}$)', size=fs)
    ax.set_ylabel('Albedo', size=fs)
    ax.set_title('{m} {v} {n}'.format(m=model, v=var_name, n=region), size=fs)
    ax.tick_params(axis='both',labelsize=fs)
    if len(olr)>10:
        ax.text(300,0.75,"{l} Profiles".format(l=len(olr)), fontsize=fs, color="0.3", ha="right")

    # plot the colorbar
    if colorbar_on:
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(csn, ax=ax, shrink=0.8)
        cb.ax.tick_params(labelsize=fs-4)
        if stat=="density":
            cb.set_label('log10(pdf)', fontsize=fs)
        elif stat=="difference":
            cb.set_label('pdf % difference', fontsize=fs)
            cb.set_ticks((levels[1:]+levels[:-1])/2)
        else:
            cb.set_label('log10(%s) (%s)'%(stat, units), fontsize=fs)
    if save:
        plt.savefig('../plots/olr_alb/jhist_%s_%s_%s_%s.png'%(var_name.lower().replace(" ","_"),
                                                               stat, model, region[:3]), bbox_inches="tight")
        print('    saved as ../plots/olr_alb/jhist_%s_%s_%s_%s.png'%(var_name.lower().replace(" ","_"),
                                                                   stat, model, region[:3]))
    if colorbar_on:
        ret = ax
    else:
        ret = ax, csn
    return ret

