# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 18:36:30 2022

@author: mugdhapolimera
"""

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
from scipy.io.idl import readsav
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
sdsscat = 'jhu'
s06file = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
bptfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
midirfile = 'mid_ir/ECO+RESOLVE_WISE_good_final.csv'
resdfname = "ECO+RESOLVE_barysample.csv"

#os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
#sdsscat = 'jhu'
#s06file = 'resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#bptfile = 'resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
#midirfile = 'mid_ir/RESOLVE_WISE_good_randerr.csv'
#resdfname = "RESOLVE_barysample.csv"


s06 = pd.read_csv(s06file)
s06.index = s06.galname
bpt = pd.read_csv(bptfile)
bpt.index = bpt.galname
midir = pd.read_csv(midirfile)
midir.index = midir.name

resdf = pd.read_csv(resdfname)
resdf.index = resdf.name
#fig, ax = plt.subplots(3,3)

bpt['agn'] = ~(bpt['defstarform'])
bpt['name'] = bpt['galname']
s06['agn'] = ~(s06['defstarform'])
s06['name'] = s06['galname']
midir['agn'] = midir['agnflag']

###############################################################################
#stdagn = np.intersect1d(bpt.galname[~(bpt['sftoagn'])], bpt.galname[bpt['agn']])
#commonstdagn = s06['agn'] & ~bpt['agn']
#s06bonusagn = np.array(s06.galname[~commonstdagn & s06['agn']])
#sfagn = np.array(bpt.galname[bpt['sftoagn']])
sf = np.array(list(bpt.galname[bpt['defstarform']]) + list(s06.galname[s06['defstarform']]))
stdagn = np.array(list(bpt.galname[bpt['composite'] | bpt['defagn'] | bpt['agntosf']]) + \
                        list(s06.galname[s06['defagn']]))
s06bonusagn = np.array(s06.galname[s06['composite']])
sfagn = np.array(bpt.galname[bpt['sftoagn']])
midiragn = np.array(midir.name[midir['agn']])


###############################################################################


keys = ['agn'] #defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
    
marker = {'agntosf': 'c^', 'ambigagn': 'ms', 'composite': 'ms', 'defagn': 'ro', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.','sftoagn': 'bs', 'sftoagn1': 's', 'sftoagn2': 'm*'}

#colors = {'agntosf': 'c', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
#          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
#          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}

labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Traditional AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'SF-AGN', 'sftoagn2' : 'MP-AGN2'}


plt.figure()

#xmin = resdf.logmstar.min(); xmax = resdf.logmstar.max()

resdf['gas_star_ratio'] = np.log10(10**resdf.logmgas/10**resdf.logmstar)
xmin = 7.5; xmax = 12.5 
ymin = 10**-3; ymax = 1.0
nbins = 25

sf_contour = np.column_stack((resdf.logmstar.loc[sf], resdf.gas_star_ratio.loc[sf]))
xmin_sf = sf_contour[:,0].min(); xmax_sf = sf_contour[:,0].max()
ymin_sf = sf_contour[:,1].min(); ymax_sf = sf_contour[:,1].max()
xgrid_sf, ygrid_sf = np.mgrid[xmin_sf:xmax_sf:nbins*1j, 
                                ymin_sf:ymax_sf:nbins*1j]
k = kde.gaussian_kde(sf_contour.T)
sf_contour_z = k(np.vstack([xgrid_sf.flatten(), ygrid_sf.flatten()]))
plt.contour(xgrid_sf, ygrid_sf, sf_contour_z.reshape(xgrid_sf.shape), 
            colors = 'k')

stdagn_contour = np.column_stack((resdf.logmstar.loc[stdagn], resdf.gas_star_ratio.loc[stdagn]))
xmin_agn = stdagn_contour[:,0].min(); xmax_agn = stdagn_contour[:,0].max()
ymin_agn = stdagn_contour[:,1].min(); ymax_agn = stdagn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(stdagn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
sf_colors_map = truncate_colormap(cm.gray_r)
plt.pcolormesh(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
                   shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
#plt.plot(resdf.logmstar.loc[stdagn],resdf.gas_star_ratio.loc[stdagn], 
#               'o', color = 'black', label = 'Standard AGN', markersize = 5,
#               alpha = 0.1, mew = 0, zorder = 0)
plt.plot(resdf.logmstar.loc[s06bonusagn],resdf.gas_star_ratio.loc[s06bonusagn], 
               '1', color = 'lime', label = 'S06 Bonus AGN', markersize = 10, mew = 2,
               alpha = 1.0, zorder = 4)
plt.plot(resdf.logmstar.loc[sfagn], resdf.gas_star_ratio.loc[sfagn], 
               's', color = 'blue', label = 'SF-AGN', markersize = 10, mew = 2, mfc = 'none',
               alpha = 1.0, zorder = 3)
plt.plot(resdf.logmstar.loc[midiragn], resdf.gas_star_ratio.loc[midiragn], 
               'p', color = 'orange', label = 'Mid-IR AGN', markersize = 10, mew = 0, zorder = 5)
plt.plot(np.linspace(7.5,12.5), 
               np.zeros(len(np.linspace(7.5,12.5))), 'k-.')
#plt.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
plt.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=15, color='k')

plt.text(9.52, -2.2, r'Gas Richness', fontsize=15, color='k')
plt.text(9.52, -2.4, r'Threshold Mass', fontsize=15, color='k')

#plt.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
#               np.linspace(10**-2.5,10**1.5), 'k-.')
plt.plot(9.5*np.ones(len(np.linspace(-2.5,2.5))), np.linspace(-2.5,2.5), 'k-.')
#plt.ylim(10**-2.5,10**1.5)
plt.ylim(-2.5,2.5)
plt.xlim(7.5,12.5) 
#plt.yscale("log")
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
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
Z_pp04[bad] = np.nan

yaxis = np.linspace(np.min(Z_pp04)-0.5,np.max(Z_pp04)+0.5)
xaxis = np.linspace(7, 11.5)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
Z04 = np.log10(0.4)+8.76
plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
#ax2 = ax1.twinx()
#ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
#float_formatter = lambda x: "%.2f" % x
#xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
#N2 = 1.754*yticks - 15.614
#N2_label = ["%.2f" % z for z in N2]
#ax2.set_yticklabels(N2_label)
#ax2.set_ylabel(r'[NII]/H$\alpha$')#, fontsize = 22)

#plt.ylim(min(yaxis),max(yaxis))
ax1.set_ylim(7.8,9.2)
plt.xlim(7.5,11.5)
ax1.set_xlabel(r'log(M${\rm_{stellar}}$/M${_\odot}$)')#, fontsize = 22)
ax1.set_ylabel(r'Z (12 + log(O/H))')#, fontsize = 22)
real = np.isfinite(Z_pp04)
realdf = resdf[real]
real_Z_pp04 = Z_pp04[real]
X,Y,Z = density_estimation(realdf.logmstar.loc[np.intersect1d(sf, realdf.index)],
                          real_Z_pp04.loc[np.intersect1d(sf, realdf.index)])
Z = Z/Z.max()

nbins = 100
#stdagn_contour = np.column_stack((realdf.logmstar.loc[np.intersect1d(stdagn, realdf.index)], 
#                                  real_Z_pp04.loc[np.intersect1d(stdagn, realdf.index)]))
stdagn_contour = np.column_stack((realdf.logmstar, real_Z_pp04))
#xmin_agn = stdagn_contour[:,0].min(); xmax_agn = stdagn_contour[:,0].max()
#ymin_agn = stdagn_contour[:,1].min(); ymax_agn = stdagn_contour[:,1].max()
xmin_agn = xmin; xmax_agn = xmax
ymin_agn = ymin; ymax_agn = ymax
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(stdagn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
agn_colors_map = truncate_colormap(cm.gray_r)
ax1.pcolormesh(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
                   shading='gouraud', cmap=agn_colors_map) #plt.cm.gray_r)

sf_colors_map = truncate_colormap(cm.gray,n=11)
nbins = 20
lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#ax1.contour(X, Y, Z, levels = lvls, colors = 'k', zorder = 1)

ax1.plot(resdf.logmstar.loc[midir.name[midir['agn']]], Z_pp04[midir.name[midir['agn']]], 
         'p', markersize = 10, label = 'Mid-IR AGN', 
         color = 'orange', alpha = 1.0, mec = 'none', zorder = 5)
ax1.plot(resdf.logmstar.loc[s06bonusagn], Z_pp04.loc[s06bonusagn], '1', markersize = 10, label = 'S06 Bonus AGN', 
         color = 'lime', alpha = 1.0, mew = 2, zorder = 4)
ax1.plot(resdf.logmstar.loc[stdagn], Z_pp04.loc[stdagn], 'o', markersize = 10, label = 'Traditional AGN', 
         color = 'blue', alpha = 1.0, mec = 'red', mew = 2, mfc = 'none', zorder = 2)
ax1.plot(resdf.logmstar.loc[sfagn], Z_pp04.loc[sfagn], 's', markersize = 10, label = 'SF-AGN', 
         color = 'blue', alpha = 1.0, mec = 'blue', mew = 2, mfc = 'none', zorder = 3)
ax1.legend(loc='lower right', fontsize = 18)


###############################################################################
#Substitute RESOLVE SFR_NUV_WISE to SEL-DF SFR_NUV from database
rescat = readsav("resolvecatalog_031622.dat",python_dict=True)
rescatphot = readsav("resolvecatalogphot_021622.dat",python_dict=True)
match, rescatndx, resselndx = np.intersect1d(rescat['name'], resdf.name, return_indices=True)
rescatsel_sfr = rescat['sfr_nuv_wise'][rescatndx]
resdf.loc[resdf.name.iloc[resselndx],'sfr_nuv_wise'] = rescatsel_sfr

ecocat = pd.read_csv("sfr_nuv_wisew_ecoSEL_newsnr.txt")
ecocat['sfr_nuv_wisew'][(ecocat['mw4w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw4w'] < 0)]
ecocat['sfr_nuv_wisew'][(ecocat['mw3w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw3w'] < 0)]
match, ecocatndx, ecoselndx = np.intersect1d(ecocat['econame'], resdf.name, return_indices=True)
ecocatsel_sfr = np.array(ecocat['sfr_nuv_wisew'][ecocatndx])
ecocatsel_sfr_err = np.array(ecocat['sfr_nuv_wisew_err'][ecocatndx])
resdf.loc[resdf.name.iloc[ecoselndx],'sfr_nuv_wise'] = ecocatsel_sfr


ssfr_st = 10**(np.log10(resdf.sfr_nuv_wise) - resdf.logmstar)
ssfr_lt = 1/(1+(1/(resdf.meanfsmgr)))
#fsmgr_st = 100*(10**6)*(ssfr_st)/(0.1*1e9*(1-ssfr_st)) #using ssfr instead of sfr
fsmgr_st = np.log10(100*(1e6)*(resdf.sfr_nuv_wise)/((10**resdf.logmstar - (100*1e6*resdf.sfr_nuv_wise)) *0.1*1e9))
fsmgr_lt = np.log10(resdf.meanfsmgr)

plt.figure()
#xaxis = np.linspace(10**-12.3, 1.2*np.max(fsmgr_st),num=50)
xaxis = np.linspace(-12.3, np.max(fsmgr_st),num=50)
yaxis = np.zeros(len(xaxis))

xmin = np.min(xaxis); xmax = np.max(xaxis)
ymin = np.min(fsmgr_lt); ymax = np.max(fsmgr_lt)
#ymin = 0.0001; ymax = 15

sf_contour = np.column_stack((fsmgr_st.loc[sf], fsmgr_lt.loc[sf]))
xmin_sf = xmin; xmax_sf = xmax
ymin_sf = ymin; ymax_sf = ymax
xgrid_sf, ygrid_sf = np.mgrid[xmin_sf:xmax_sf:nbins*1j, 
                                ymin_sf:ymax_sf:nbins*1j]
k = kde.gaussian_kde(sf_contour.T)
sf_contour_z = k(np.vstack([xgrid_sf.flatten(), ygrid_sf.flatten()]))
nbins = 20
lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
plt.contour(xgrid_sf, ygrid_sf, sf_contour_z.reshape(xgrid_sf.shape), 
            colors = 'k')

nbins = 30
stdagn_contour = np.column_stack((fsmgr_st.loc[stdagn], fsmgr_lt.loc[stdagn]))
xmin_agn = xmin; xmax_agn = xmax
ymin_agn = ymin; ymax_agn = ymax
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(stdagn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
agn_colors_map = truncate_colormap(cm.gray_r, n=5)
plt.pcolormesh(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
                   shading='gouraud', cmap=agn_colors_map) 

plt.plot(fsmgr_st.loc[s06bonusagn], fsmgr_lt.loc[s06bonusagn], 
               '1', color = 'lime', label = 'S06 AGN', markersize = 10, mew = 2,
               alpha = 1, zorder = 4)
plt.plot(fsmgr_st.loc[sfagn], fsmgr_lt.loc[sfagn], 's', color = 'blue', 
         mec = 'blue', mfc = 'none', label = 'Optimized scheme SF-AGN', 
         markersize = 10, mew = 2, alpha = 1, zorder = 5)
plt.plot(fsmgr_st.loc[midiragn], fsmgr_lt.loc[midiragn], 
               'p', color = 'orange', label = 'Mid-IR AGN', markersize = 10, mew = 0, zorder = 5)
plt.plot(xaxis, yaxis, 'k--', lw = 3)
plt.xlim(np.min(xaxis), np.max(fsmgr_st))
plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
plt.text(0.0005, 1.25, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
plt.tick_params(which = 'minor', length=4, width=2)

#plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
#            '*',color = 'lime', ms = 10, label = labels[key])
plt.xlabel(r'Short Term SFH $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr) \times 0.1 Gyr}\right)$')
plt.ylabel(r'Long Term SFH $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr) \times 1 Gyr}\right)$')

###############################################################################

#plt.figure()
#sel = resdf.loc[s06.name[s06['agn']]]
#mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
#plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], '1', color = 'red',
#         markersize = 5, alpha = 1.0, mew = 2) #other group galaxies
#plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], '1', color = 'red',
#         markersize = 10, alpha = 1.0, zorder = 2, mew = 2, label = 'S06 AGN') #central
#sel = resdf.loc[bpt.name[bpt['sftoagn']]]
#mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
#plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 's', color = 'blue',
#         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
#plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 's', color = 'blue',
#         markersize = 10, alpha = 1.0, mec = 'none', zorder = 3, label = 'Optimized scheme SF-AGN') #central
#sel = resdf.loc[bpt.name[bpt['agn']]]
#mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
#plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 's', color = 'green',
#         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
#plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 's', color = 'green',
#         markersize = 10, alpha = 1.0, mec = 'none', zorder = 0, label = 'Optimized scheme other AGN') #central
#sel = resdf.loc[midir.name[midir['agn']]]
#mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
#plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 'p', color = 'orange',
#         markersize = 5, alpha = 1.0, mec = 'none') #other group galaxies
#plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 'p', color = 'orange',
#         markersize = 10, alpha = 1.0, mec = 'none', zorder = 4, label = 'Mid-IR AGN') #central
#
#plt.legend(loc = 'lower right', numpoints = 1, fontsize = 15)
#plt.xlabel(r'log(Galaxy M$\rm_{halo}$/M$_{\odot}$)')
#plt.ylabel(r'log(Galaxy M$\rm_{bary}$/M$_{\odot}$)')
#plt.xlim(10.8,14.54)
#plt.ylim(9.15,11.5)
#fullmbary = np.log10(10**df.logmstar + 10**df.logmgas)
#yaxis = np.arange(np.min(fullmbary)-0.5, np.max(fullmbary)+0.5,0.1)
#xaxis = np.ones(len(yaxis))
#plt.plot(11.5*xaxis, yaxis, 'k--', lw = 3)
#plt.plot(12.0*xaxis, yaxis, 'k-.', lw = 3)

for x in resdf.keys():
    resdf[x] = np.array(resdf[x]).byteswap().newbyteorder() 

xmin = 7.5
xmax = 11.5
ymin = 0
ymax = 3
u_r = resdf['modelu_rcorr']

X,Y,Z = density_estimation(resdf.logmstar,u_r)
fig,ax = plt.subplots()#(figsize=(8,8))
ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
sel = resdf.loc[np.array(s06[s06['agn']].name)]
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
