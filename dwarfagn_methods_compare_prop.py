# -*- coding: utf-8 -*-
"""
Comparison of properties probed by different AGN detection methods
Created on Thu Jun 17 15:57:27 2021

@author: mugdhapolimera
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colors as mpcolors
import scipy
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
sdsscat = 'jhu'
s06file = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
bptfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
midirfile = 'mid_ir/RESOLVE_WISE_good.csv'

s06 = pd.read_csv(s06file)
s06.index = s06.galname
bpt = pd.read_csv(bptfile)
bpt.index = bpt.galname
midir = pd.read_csv(midirfile)
midir.index = midir.name

resdfname = "ECO+RESOLVE_barysample.csv"
resdf = pd.read_csv(resdfname)
resdf.index = resdf.name
#fig, ax = plt.subplots(3,3)

bpt['agn'] = ~(bpt['defstarform'])
bpt['name'] = bpt['galname']
s06['agn'] = ~(s06['defstarform'])
s06['name'] = s06['galname']
midir['agn'] = midir['agnflag']

keys = ['agn'] #defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
    
marker = {'agntosf': 'c^', 'ambigagn': 'ms', 'composite': 'ms', 'defagn': 'ro', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.','sftoagn': 'bs', 'sftoagn1': 's', 'sftoagn2': 'm*'}

colors = {'agntosf': 'c', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}

labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Traditional AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'SF-AGN', 'sftoagn2' : 'MP-AGN2'}


plt.figure()
plt.plot(resdf.logmstar.loc[s06.name[s06['agn']]],
               10**resdf.logmgas.loc[s06.name[s06['agn']]]/10**resdf.logmstar.loc[s06.name[s06['agn']]], 
               '1', color = 'red', label = 'S06 AGN', markersize = 10, mew = 2,
               alpha = 1.0, zorder = 4)
plt.plot(resdf.logmstar.loc[bpt.name[bpt['sftoagn']]],
               10**resdf.logmgas.loc[bpt.name[bpt['sftoagn']]]/10**resdf.logmstar.loc[bpt.name[bpt['sftoagn']]], 
               's', color = 'blue', label = 'Optimized scheme SF-AGN', markersize = 10, mew = 0,
               alpha = 1.0, zorder = 3)
plt.plot(resdf.logmstar.loc[bpt.name[bpt['agn']]],
               10**resdf.logmgas.loc[bpt.name[bpt['agn']]]/10**resdf.logmstar.loc[bpt.name[bpt['agn']]], 
               's', color = 'green', label = 'Optimized scheme Other AGN', markersize = 10, mew = 0,
               alpha = 1.0)
plt.plot(resdf.logmstar.loc[midir.name[midir['agn']]],
               10**resdf.logmgas.loc[midir.name[midir['agn']]]/10**resdf.logmstar.loc[midir.name[midir['agn']]], 
               'p', color = 'orange', label = 'Mid-IR AGN', markersize = 10, mew = 0, zorder = 5)
plt.plot(np.linspace(7.5,12.5), 
               np.ones(len(np.linspace(7.5,12.5))), 'k-.')
#plt.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
plt.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=15, color='k')

plt.text(9.52, -2.2, r'Gas Richness', fontsize=15, color='k')
plt.text(9.52, -2.4, r'Threshold Mass', fontsize=15, color='k')

plt.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
               np.linspace(10**-2.5,10**1.5), 'k-.')
plt.ylim(10**-2.5,10**1.5)
plt.xlim(7.5,12.5) 
plt.yscale("log")
plt.legend(loc='upper right',numpoints = 1, fontsize = 15) # bbox_to_anchor=(1,1),
plt.xlabel(r'$\rm \log(M_{stellar}/M_{\odot})$')
plt.ylabel(r'$\rm M_{gas}/M_{stellar}$')

###############################################################################
plt.figure()
ymin, ymax = (7.8, 9.2)
xmin,xmax = (7.5,11.5)
ax1 = plt.subplot()
N2 = np.log10(resdf.nii_6584_flux/resdf.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
df = resdf[(~np.isnan(N2))]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
yaxis = np.linspace(np.min(Z_pp04)-0.5,np.max(Z_pp04)+0.5)
xaxis = np.linspace(7, 11.5)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
Z04 = np.log10(0.4)+8.76
plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
ax2 = ax1.twinx()
yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
ax2.set_yticks(np.linspace(0,1,len(yticks)))
#ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
float_formatter = lambda x: "%.2f" % x
#xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
N2 = 1.754*yticks - 15.614
N2_label = ["%.2f" % z for z in N2]
ax2.set_yticklabels(N2_label)
ax2.set_ylabel(r'[NII]/H$\alpha$')#, fontsize = 22)

#plt.ylim(min(yaxis),max(yaxis))
ax1.set_ylim(7.8,9.2)
plt.xlim(7.5,11.5)
ax1.set_xlabel(r'log(M${\rm_{stellar}}$/M${_\odot}$)')#, fontsize = 22)
ax1.set_ylabel(r'Z (12 + log(O/H))')#, fontsize = 22)
X,Y,Z = density_estimation(df.logmstar,Z_pp04)
Z = Z/Z.max()
#    ax1.imshow(np.rot90(Z), cmap='bone_r',                                                    
#              extent=[xmin, xmax, ymin, ymax], interpolation='gaussian')
#    plt.clim(0,1.8)
lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#    lvls = [ 0.25, 0.5, 0.75, 1.0]
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = mpcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
sf_colors_map = truncate_colormap(cm.gray,n=11)
nbins = 20
ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 4)



df = resdf.loc[midir.name]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
df = df[(~np.isnan(N2))]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
ax1.plot(df.logmstar[midir.agn], Z_pp04[midir.agn], 'p', markersize = 10, label = 'Mid-IR AGN', 
         color = 'orange', alpha = 1.0, mec = 'none', zorder = 5)
df = resdf.loc[s06.name]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
df = df[(~np.isnan(N2))]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
ax1.plot(df.logmstar[s06.agn], Z_pp04[s06.agn], '1', markersize = 10, label = 'S06 AGN', 
         color = 'red', alpha = 1.0, mew = 2, zorder = 4)
df = resdf.loc[bpt.name[bpt.sftoagn]]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
df = df[(~np.isnan(N2))]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
ax1.plot(df.logmstar[bpt.sftoagn], Z_pp04[bpt.sftoagn], 's', markersize = 10, label = 'Optimized scheme SF-AGN', 
         color = 'blue', alpha = 1.0, mec = 'none', zorder = 3)
df = resdf.loc[bpt.name]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
df = df[(~np.isnan(N2))]
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
ax1.plot(df.logmstar[bpt.agn], Z_pp04[bpt.agn], 's', markersize = 10, label = 'Optimized scheme Other AGN', 
         color = 'green', alpha = 1.0, mec = 'none')
ax1.legend(loc='lower right', fontsize = 18)


###############################################################################
ssfr_st = 10**(np.log10(resdf.sfr_nuv) - resdf.logmstar)
ssfr_lt = 1/(1+(1/(resdf.meanfsmgr)))
#fsmgr_st = 100*(10**6)*(ssfr_st)/(0.1*1e9*(1-ssfr_st)) #using ssfr instead of sfr
fsmgr_st = 100*(1e6)*(resdf.sfr_nuv)/((10**resdf.logmstar - (100*1e6*resdf.sfr_nuv)) *0.1*1e9)
fsmgr_lt = resdf.meanfsmgr

plt.figure()
plt.plot(fsmgr_st.loc[s06.name[s06['agn']]],
               fsmgr_lt.loc[s06.name[s06['agn']]], 
               '1', color = 'red', label = 'S06 AGN', markersize = 10, mew = 2,
               alpha = 1, zorder = 4)
plt.plot(fsmgr_st.loc[bpt.name[bpt['sftoagn']]],
               fsmgr_lt.loc[bpt.name[bpt['sftoagn']]], 
               's', color = 'blue', label = 'Optimized scheme SF-AGN', markersize = 10, mew = 0,
               alpha = 1, zorder = 3)
plt.plot(fsmgr_st.loc[bpt.name[bpt['agn']]],
               fsmgr_lt.loc[bpt.name[bpt['agn']]], 
               's', color = 'green', label = 'Optimized scheme Other AGN', markersize = 10, mew = 0,
               alpha = 1.0)
plt.plot(fsmgr_st.loc[midir.name[midir['agn']]],
               fsmgr_lt.loc[midir.name[midir['agn']]], 
               'p', color = 'orange', label = 'Mid-IR AGN', markersize = 10, mew = 0, zorder = 5)
plt.plot(np.linspace(7.5,12.5), 
               np.ones(len(np.linspace(7.5,12.5))), 'k-.')
xaxis = np.arange(10**-12.3, np.max(fsmgr_st)+0.05,0.01)
yaxis = np.ones(len(xaxis))
plt.plot(xaxis, yaxis, 'k--', lw = 3)
plt.xlim(np.min(xaxis), np.max(fsmgr_st))
plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
plt.text(0.0005, 1.25, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
plt.tick_params(which = 'minor', length=4, width=2)

#plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
#            '*',color = 'lime', ms = 10, label = labels[key])
plt.xlabel(r'Short Term FSMGR $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr) \times 0.1 Gyr}\right)$')
plt.ylabel(r'Long Term FSMGR $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr) \times 1 Gyr}\right)$')

###############################################################################

plt.figure()
sel = resdf.loc[s06.name[s06['agn']]]
mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], '1', color = 'red',
         markersize = 5, alpha = 1.0, mew = 2) #other group galaxies
plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], '1', color = 'red',
         markersize = 10, alpha = 1.0, zorder = 2, mew = 2, label = 'S06 AGN') #central
sel = resdf.loc[bpt.name[bpt['sftoagn']]]
mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 's', color = 'blue',
         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 's', color = 'blue',
         markersize = 10, alpha = 1.0, mec = 'none', zorder = 3, label = 'Optimized scheme SF-AGN') #central
sel = resdf.loc[bpt.name[bpt['agn']]]
mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 's', color = 'green',
         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 's', color = 'green',
         markersize = 10, alpha = 1.0, mec = 'none', zorder = 0, label = 'Optimized scheme other AGN') #central
sel = resdf.loc[midir.name[midir['agn']]]
mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 'p', color = 'orange',
         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 'p', color = 'orange',
         markersize = 10, alpha = 1.0, mec = 'none', zorder = 4, label = 'Mid-IR AGN') #central

plt.legend(loc = 'lower right', numpoints = 1, fontsize = 15)
plt.xlabel(r'log(Galaxy M$\rm_{halo}$/M$_{\odot}$)')
plt.ylabel(r'log(Galaxy M$\rm_{bary}$/M$_{\odot}$)')
plt.xlim(10.8,14.54)
plt.ylim(9.15,11.5)
fullmbary = np.log10(10**df.logmstar + 10**df.logmgas)
yaxis = np.arange(np.min(fullmbary)-0.5, np.max(fullmbary)+0.5,0.1)
xaxis = np.ones(len(yaxis))
plt.plot(11.5*xaxis, yaxis, 'k--', lw = 3)
plt.plot(12.0*xaxis, yaxis, 'k-.', lw = 3)



xmin = 7.5
xmax = 11.5
ymin = 0
ymax = 3
u_r = resdf['modelu_rcorr']

X,Y,Z = density_estimation(resdf.logmstar,u_r)
fig,ax = plt.subplots()#(figsize=(8,8))
ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
sel = resdf.loc[s06.name[s06['agn']]]
ax.plot(sel.logmstar, sel.modelu_rcorr, '1', color = 'red', markersize = 10, 
       mew = 2, label = 'S06 AGN', zorder = 3)
sel = resdf.loc[bpt.name[bpt['sftoagn']]]
ax.plot(sel.logmstar, sel.modelu_rcorr, 's', color = 'blue', markersize = 10, 
        label = 'Optimized scheme SF-AGN', zorder = 2)
sel = resdf.loc[bpt.name[bpt['agn']]]
ax.plot(sel.logmstar, sel.modelu_rcorr, 's', color = 'green', markersize = 10, 
        label = 'Optimized scheme other AGN', zorder = 1)
sel = resdf.loc[midir.name[midir['agn']]]
ax.plot(sel.logmstar, sel.modelu_rcorr, 'p', color = 'orange', markersize = 10, 
        label = 'Mid-IR AGN', zorder = 5)
ax.axhline(y = 1.4)
#ax.legend(handles=[sfagnplt,targetsplt,samiplt], 
#          labels = ['SFing-AGN', 'Targets for this proposal','SFing-AGN in SAMI'],
#          fontsize = 12, loc = 'upper left')
plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
plt.ylabel('u-r', fontsize = 15)
ax.contour(X, Y, Z, cmap='summer')
#
