# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab
import os
import DP_diagnostics_functions

def plot1Dtime(datadir,plotdir,time_start,time_end,filelist,caselist,\
               xaxis_opt="default",xaxis_units="time (days)",\
               xaxis_mult=0,xaxis_start=0):
    
    colorarr=["r","b","g","c","m","y"]
    
    numfiles=len(filelist)
    
    ############################################################
    # Make the variable list for variables to plot
    # This will make profile plots for every variable possible
    
    # Initialize list of variables to plot
    varstoplot=[];
    
    # temporarily limit to searching just the first file for speed
    listtodo=filelist[0]
    DP_diagnostics_functions.makevarlist(datadir,listtodo,3,varstoplot)
    
    print('Now Plotting 1D Time Fields')
    ############################################################
    # loop over the variables to plot
    for x in range(0,len(varstoplot)):
        
        varname=varstoplot[x]
        legendlist=[];
    
        # loop over the number of simulations to plot
        for f in range(0,numfiles):
            
            filename=filelist[f]
            file=datadir+filename
            
            # Read in file and the dataset information
            fh=Dataset(file,mode='r')
            
            # Read in appropriate coordinate variables
            time=fh.variables['time'][:]
            
            # Generate a list to see if variables are in there
            varsinfile=fh.variables.keys()
                
            if varname in varsinfile:
                vartoplot=fh.variables[varname][:]
                
                if (time_end == "end"):
                    plottimes=np.squeeze(np.where(time >= time_start))
                else:
                    plottimes=np.squeeze(np.where((time >= time_start) & \
                                                  (time <= time_end)))

                if (vartoplot.ndim == 3):                    
                    vartoplot2=np.squeeze(vartoplot[plottimes,:,:])
                elif (vartoplot.ndim == 2):
                    vartoplot2=np.mean(vartoplot[plottimes,:],axis=1)

                if (xaxis_opt == "custom"):
                    time=(time*xaxis_mult)+xaxis_start
                else:
                    xaxis_units="time (days)"
            
                plt.figure(x)
                plt.plot(time[plottimes],vartoplot2,colorarr[f],linewidth=3)
                
                legendlist.append(caselist[f])
        
                plt.title(fh.variables[varname].long_name + '\n' + \
                '(' + varname + ')',fontsize=16)
                plt.xlabel(xaxis_units,fontsize=14)
                if hasattr(fh.variables[varname],'units'):
                    plt.ylabel('('+fh.variables[varname].units+')',fontsize=14)
                plt.grid(True)
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                
        plt.legend(legendlist)       
        pylab.savefig(plotdir+'/'+varname+'.png',bbox_inches='tight',pad_inches=0.05)
        plt.close(x)

    return varstoplot

################################################################
### Plot the 2 dimensional timeseries

def plot2Dtime(datadir,plotdir,toplev,filelist,caselist): 
    
    numfiles=len(filelist)
    
    # Make list of variables to plot
    varstoplot=[];
    DP_diagnostics_functions.makevarlist(datadir,filelist,4,varstoplot)

    print('Now Plotting 2D Time Fields')
    ############################################################
    # loop over the variables to plot
    for x in range(0,len(varstoplot)):
        
        varname=varstoplot[x]
    
        # loop over the number of simulations to plot
        for f in range(0,numfiles):
            
            filename=filelist[f]
            file=datadir+filename
            
            # Read in file and the dataset information
            fh=Dataset(file,mode='r')
            
            # Read in appropriate coordinate variables
            time=fh.variables['time'][:]
            lev=fh.variables['lev'][:]
            ilev=fh.variables['ilev'][:]
            
            # Generate a list to see if variables are in there
            varsinfile=fh.variables.keys()
                
            if varname in varsinfile:
                vartoplot=fh.variables[varname][:]
                
                avgprof=np.mean(vartoplot,axis=0)
                
                if (len(avgprof) == len(lev)):
                    levarr=lev
                elif (len(avgprof) == len(ilev)):
                    levarr=ilev
            
                # Only plot to top level
                plotlevs=np.where(levarr > toplev)

                plt.figure(x)
                vartoplot=np.squeeze(vartoplot[:,plotlevs])
                vartoplot=np.rot90(vartoplot)
                levarr=levarr[plotlevs]
                rev_lev=levarr[::-1]
                sc=plt.contourf(time,rev_lev,vartoplot)
                plt.ylim(max(rev_lev),toplev)
        
                plt.title(fh.variables[varname].long_name + \
                ' (' + varname + ')' '\n' + caselist[f])
                plt.ylabel('P (hPa)')
                plt.xlabel('time (days)')
                clb=plt.colorbar(sc)
                clb.set_label('('+fh.variables[varname].units+')')
                
                pylab.savefig(plotdir+'/'+varname+caselist[f]+'.png')
                plt.close(x)
                
    return

################################################################
### Plot the 2 dimensional timeseries

def plot2Dtimepanel(datadir,plotdir,toplev,filelist,caselist,\
               xaxis_opt="default",xaxis_units="time (days)",\
               xaxis_mult=0,xaxis_start=0):
    
    numfiles=len(filelist)
    
    # Make list of variables to plot
    varstoplot=[];
    listtodo=filelist[0];
    DP_diagnostics_functions.makevarlist(datadir,listtodo,4,varstoplot)

#    varstoplot=["CLOUD"]

    print('Now Plotting 2D Time Fields')
    ############################################################
    # loop over the variables to plot
    for x in range(0,len(varstoplot)):
        
        varname=varstoplot[x]
        varsinpanel=0
    
        # loop over the number of simulations to plot
        for f in range(0,numfiles):
            
            filename=filelist[f]
            file=datadir+filename
            
            # Read in file and the dataset information
            fh=Dataset(file,mode='r')
            
            # Generate a list to see if variables are in there
            varsinfile=fh.variables.keys()
                
            if varname in varsinfile:      
                varsinpanel=varsinpanel+1

        # start loop again and make in panel
         # loop over the number of simulations to plot
         
        plt.figure(x)
        if varsinpanel > 0: 
            if varsinpanel == 1:
                fig,axes=plt.subplots()
            elif (varsinpanel > 1) & (varsinpanel <= 4):
                fig,axes=plt.subplots(varsinpanel,sharex=True)
            elif (varsinpanel > 4) & (varsinpanel % 2 == 0):
                fig,axes=plt.subplots(varsinpanel/2.,2,sharey=True,sharex=True)
            elif (varsinpanel > 4) & (varsinpanel % 2 != 0):
                fig,axes=plt.subplots(varsinpanel/2.+1,2,sharey=True,sharex=True)
                
#            fig.tight_layout()
            
            pcount=0
#            plt.figure(x) 
            for f in range(0,numfiles):

                Rd=287.15 # gas constant for air
                grav=9.81 # gravity
                
                filename=filelist[f]
                file=datadir+filename
                
                # Read in file and the dataset information
                fh=Dataset(file,mode='r')
                
                # Read in appropriate coordinate variables
                time=fh.variables['time'][:]
                lev=fh.variables['lev'][:]
                ilev=fh.variables['ilev'][:]

                ###########################################
	        # Compute heights
                temp=fh.variables['T'][:]

                z_mid=np.zeros(lev.size)
                z_int=np.zeros(ilev.size)

                T_avg=np.mean(temp,axis=2)

	        # Guess the bottom edge
                z_int[-1]=0.0
                z_mid[-1]=10.0
                for k in range(len(lev)-2,-1,-1):
                    z_mid[k] = z_mid[k+1] + Rd*T_avg[0,k]/grav*np.log(lev[k+1]/lev[k])

                for k in range(len(ilev)-2,-1,-1):
                    z_int[k] = z_int[k+1] + Rd*T_avg[0,k]/grav*np.log(ilev[k+1]/ilev[k])

                # Convert from meters to kilometers
                z_mid=z_mid/1000.
                z_int=z_int/1000.

                # End compute heights
		###########################################3
                
                # Generate a list to see if variables are in there
                varsinfile=fh.variables.keys()
                    
                if varname in varsinfile:
                    vartoplot=fh.variables[varname][:]  
                    
                    avgprof=np.mean(vartoplot,axis=0)
                    # Horizontal average
                    vartoplot=np.mean(vartoplot,axis=2)
                    
                    if (len(avgprof) == len(lev)):
                        levarr=z_mid
                    elif (len(avgprof) == len(ilev)):
                        levarr=z_int
                
                    # Only plot to top level
                    plotlevs=np.where(levarr < toplev)
                       
                    vartoplot=np.squeeze(vartoplot[:,plotlevs])
                    vartoplot=np.rot90(vartoplot)
                    vartoplot=np.flipud(vartoplot)
                    levarr=levarr[plotlevs]
                    rev_lev=levarr

                    if (xaxis_opt == "custom"):
                       time=(time*xaxis_mult)+xaxis_start
                    else:
                       xaxis_units="time (days)"

                    if pcount == 0:
                        fig.suptitle(fh.variables[varname].long_name + \
                    ' (' + varname + ')',y=1.01)
                    if varsinpanel == 1:
                        sc=axes.contourf(time,rev_lev,vartoplot)
                        axes.set_ylim(0,toplev)
                        axes.set_title(caselist[f])
                        axes.set_ylabel('height (km)')
                    elif (varsinpanel > 1) & (varsinpanel <= 4):
                        sc=axes[pcount].contourf(time,rev_lev,vartoplot)
                        axes[pcount].set_ylim(0,toplev)
                        axes[pcount].set_title(caselist[f])
                        axes[pcount].set_ylabel('height (km)')
                    elif (varsinpanel > 4) & (pcount % 2 == 0):
                        sc=axes[pcount/2,0].plt.contourf(time,rev_lev,vartoplot)
                        axes[pcount/2,0].set_ylim(0,toplev)
                        axes[pcount/2,0].set_title(caselist[f])
                        axes[pcount/2,0].set_ylabel('height (km)')
                    elif (varsinpanel > 4) & (pcount % 2 != 0):
                        sc=axes[pcount/2,1].plt.contourf(time,rev_lev,vartoplot)
                        axes[pcount/2,1].set_ylim(0,toplev)
                        axes[pcount/2,1].set_title(caselist[f])
                        axes[pcount/2,1].set_ylabel('height (km)')
                        
                    # get unit info here because next file may not have it
                    if hasattr(fh.variables[varname],'units'):
                        theunits=fh.variables[varname].units
                    else:
                        theunits="-"
                    
                    pcount=pcount+1
            
#            plt.ylabel('P (hPa)')
            plt.xlabel(xaxis_units)
            
            fig.tight_layout()

            if (varsinpanel > 1):
                clb=fig.colorbar(sc,ax=axes.ravel().tolist())
            else:
                clb=fig.colorbar(sc)

            clb.set_label(theunits)
            
#            fig.subplots_adjust(top=0.9,left=0.1,bottom=0.1)
            
            pylab.savefig(plotdir+'/'+varname+'.png',bbox_inches='tight')
            plt.close(x)                    
                
    return varstoplot
