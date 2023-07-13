import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from matplotlib import colors,ticker
import matplotlib.gridspec as gridspec

#make path to the module active
import sys
sys.path.append('/glade/u/home/blazg/scripts/modules_python/')
import sam #for module ncopen
import util


path = '/glade/campaign/univ/uwas0108/blazg/SAM/OUT_3D/'
timeind= slice(None)
levind = slice(None)
lonind = slice(None)
latind = slice(None)

varnames = ['time','y','x','z']

simul = ['2cat_d3'] #1Cat
simul = ['frzproc']
print(simul)
cat = 1
cpu =''
#ybins = np.arange(173,300, 4) #-90 to 0 deg C
ybins = 10**np.linspace(-5,2,71) 

varloop = ['rei']#,'omega','convf']

hist_var = {}
hist_freq = {}
hist_counts={}
allcounts = np.zeros((70,70))

for c,car in enumerate(varloop):
    
    if c==0:
        xbins = np.arange(0,141,2)
        
    #hist_var[car] = np.zeros((len(xbins)-1,len(ybins)-1))
    #hist_freq[car] = np.zeros((len(xbins)-1,len(ybins)-1))
    hist_counts = np.zeros((3,len(xbins)-1,len(ybins)-1))

################################################################################
# ###############################################################################
varp=['p']
#time loop starts here
for s, sim in enumerate(simul): #fake loop!!!!!
    print(s,sim)

    if cat==1:
        alpha = ['A']
    if cat==2:
        alpha = ['A','B']
    if cat==6:
        cpu = '_432cpu'
        alpha = ['A','B','C','D','E','F']
        #cpu =''
    varnames3d = ['QP','QN','QV','TABS','QR','QC','NC','NR']
    varnames3d = varnames3d+ ['QI'+alpha[i] for i in range(len(alpha))]+ ['NI'+alpha[i] for i in range(len(alpha))] # +\
                   # ['DGEICE'+alpha[i]  for i in range(len(alpha))]

    data_3d = sam.ncopen(path+'walker_1cat_3840x192_1x1km_dc_frzlim3x_192.bin2D_1.nc',varnames3d,ind=timeind)
    datap = sam.ncopen(path+'walker_1cat_3840x192_1x1km_dc_frzlim3x_192.bin2D_1.nc',varp) #pressure always the same, I think

    #data_3d = sam.ncopen(path+'walker_3d_'+simul[0]+'_3888x36x192_200x200m'+cpu+'_dc_432_00000'+step+'.nc',varnames3d,ind=timeind)
    #datap = sam.ncopen(path+'walker_3d_'+simul[0]+'_3888x36x192_200x200m'+cpu+'_dc_432_0000028800.nc',varp) #pressure always the same, I think

    #if s==0: #add also QC and QR => total water, as the obs don't distinguish!
    iwc=np.nansum([data_3d['QI'+alpha[i]] for i in range(len(alpha))],axis=0) + data_3d['QC'].squeeze() + data_3d['QR'].squeeze()
    icnc=np.nansum([data_3d['NI'+alpha[i]] for i in range(len(alpha))],axis=0) +data_3d['NC'].squeeze() + data_3d['NR'].squeeze()
    iwc = np.squeeze(iwc)
    icnc= np.squeeze(icnc)
    pres0=datap['p']
    vap = data_3d['QV'].squeeze()
    
    temp=data_3d['TABS'].squeeze()
    
    pres = np.tile(pres0[np.newaxis,1,np.newaxis],(temp.shape[0],1,temp.shape[2]))
    rho = 100*pres/(287*temp*(1+1e-3*vap/0.622)/(1+1e-3*vap))

    rei = (3*(iwc*rho)/(4*0.92*np.pi*rho*icnc) )**(1/3) * 100 #see Kramer et al., 2020, page 1: mean mass radius
    rei = np.squeeze(rei)
    rei[iwc<1e-5]  = np.nan
    icnc[iwc<1e-5] = np.nan
    #iwc[iwc<1e-5]  = np.nan
    
    rei[temp>233.15] = np.nan
    icnc[temp>233.15] = np.nan

    print(np.nanmedian(rei),'rei median')
    ################################
    #data_3d = (rei)

    icnc = icnc*rho
    
    for c,car in enumerate(varloop):

        print(car)
        plotvar_x = rei.ravel()
        plotvar_y = icnc.ravel()

        if c==0:
            xbins = np.arange(0,141,2)
            xplot4=(xbins[:-1]+xbins[1:])/2.
            
        #hist_freq[car]=np.histogram2d(plotvar_x, plotvar_y, bins=[xbins,ybins], normed=False)[0] 
        hist_counts[c] = np.histogram2d(plotvar_x, plotvar_y, bins=[xbins,ybins], normed=False)[0] 
        #hist_norm1[car]= hist_norm1[car]+hist_freq[car] #sum over timesteps
        #hist_norm1= hist_norm1 + hist_freq[car] 
    allcounts = allcounts +    hist_counts
write=1
######################################################################################
# #####################SAVE VARIABLES in an NC FILE####################################
yplot =(ybins[:-1]+ybins[1:])/2.

      
if (write == 1):
    ##############WRITE IN A NC FILE##################################################

    ############################################################################
    from netCDF4 import Dataset as NetCDFFile
    
    savepath ='/glade/campaign/univ/uwas0108/blazg/SAM/processed/icnc_rei_space/'
    
    ncfile = NetCDFFile(savepath+'ICNCvsREI_Walker2D_1km_qi1e-5_frzlim3x_Tless-40.nc', 'w',  format='NETCDF4_CLASSIC') 
    
    icnc = ncfile.createDimension('icnc', len(yplot) )
    rad = ncfile.createDimension('rad', len(xplot4) )

    # Create coordinate variables for 2-dimensions
    rad = ncfile.createVariable('ice radius', np.float32, ('rad',)) 
    icnc = ncfile.createVariable('ice number concentration', np.float32, ('icnc',)) 

    # Create the actual 2-d variable
    rei_bin    = ncfile.createVariable('binned in-cloud ice crystal radius', np.dtype('float32').char, ('icnc','rad'))

    # Variable Attributes  
    rei_bin.units      = 'micron'  
   
    #assign values
    icnc[:]           = yplot
    rad[:]            = xplot4

    rei_bin[:]       = np.flipud(allcounts[0][:,::-1].T)

    ncfile.close()        
