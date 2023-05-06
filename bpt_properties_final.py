# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:25:30 2018

@author: mugdhapolimera

This code explores the properties of galaxies categorized using BPT plots.
"""

import numpy as np
import os
import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter
from scipy.io.idl import readsav
label_size = 15
import scipy
import matplotlib as mpl
mpl.rcParams.update({'font.size': 20})
mpl.rcParams.update({'axes.linewidth': 2})
mpl.rcParams.update({'lines.linewidth': 2})
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.width'] = 2
#from scipy.stats import norm

#path = os.path.dirname(os.path.realpath(__file__))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
he2_flag = 0
full = 1
resolve = 0
eco = 0
sdsscat = 'jhu'#'nsa'
if he2_flag:
    flags = pd.read_csv('eco+resolve_emlineclass_filter_he2.csv')
else:
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if full : 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_jhu.csv'
    flags = pd.read_csv(r'ECO\SEL\eco+resolve_emlineclass_dext_snr5_jhu.csv')
if resolve: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_dext_'+sdsscat+'.csv'
    flags = pd.read_csv('resolve_emlineclass_dext_snr5_'+sdsscat+'.csv')
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO/SEL/ECO_full_snr5_dext_jhu.csv'
    flags = pd.read_csv(r'ECO\SEL\eco_emlineclass_dext_snr5_jhu.csv')

full_df = pd.read_csv(inputfile)
full_df.index = full_df.name
df = full_df.loc[list(flags['galname'])]
if 'NAME' not in df.keys():
    df['NAME'] = df['name']

keys = ['defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']

if he2_flag:
    keys.append('heiisel')

    
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

percent = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 57,
          'defstarform': 0, 'sftoagn': 0}

for key in keys:
    percent[key] = len(np.where(flags[key])[0])

df = df[df.logmstar > 0]
flags.index = flags.galname

jwst = ['rs1126', 'rf0053', 'rs0463', 'rf0399', 
        'rs0010', 'rs0472', 'rf0503', 'rs1195',
        'rf0013', 'rf0372', 'rs0462', 'rs1214',
        'rf0002','rs0582',
        'rs0029','rs1116']

jwst = []
#FSMGR vs. G/S
#plt.figure()
#plt.suptitle(sdsscat+' catalog')
#jwst = ['rs0010']
#for key in keys:
#    sel = df.iloc[np.where(flags[key])[0]]
#    if key == 'defstarform':
#        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
#             marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
#             color = colors[key], label = labels[key])
#    elif key == 'agntosf': 
#        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
#             marker[key], markersize = 10, mew = 1, color = colors[key],
#             mec = 'y', label = labels[key])
#    else:
#        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
#             marker[key], markersize = 10, mew = 0, color = colors[key],
#             label = labels[key])
#    plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
#    plt.plot(np.linspace(-2.1,1.6), 1+0*np.linspace(-2.1,1.6), 'k-.')
#    plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
#             rotation = 'vertical')
#    plt.text(-2.0, 1.5, r'Stellar Mass Doubled in last Gyr', 
#             fontsize=14, color='k')
#    
#    plt.xlabel(r'$\rm \log (M_{gas}/M_{stellar})$', size = 22)
#    plt.ylabel('Mean FSMGR', size = 22)
#    plt.yscale('log')
#    yticks = plt.yticks()[0]
#    plt.yticks(yticks, np.around(yticks,2))
#    plt.ylim(10**-3, 10**2)
#    plt.xlim(-2.1,1.6)
#    #plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)
#    

#G/S vs. M_star with histogram
'''rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHisty = plt.axes(rect_histy, xscale = 'log')
nullfmt = NullFormatter() # no labels
axHisty.yaxis.set_major_formatter(nullfmt)
'''
plt.figure('Physical Properties')
#plt.suptitle(sdsscat+' catalog')
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx, yscale = 'log')
nullfmt = NullFormatter() # no labels
axHistx.xaxis.set_major_formatter(nullfmt)

# the scatter plot:
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], label = labels[key], markersize = 8, 
                           alpha = 0.3, mew = 0)
        
        #elif key == 'sftoagn':
        #    axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
        #                   label = labels[key], markersize = 15, mew = 0)
                
        elif key == 'heiisel':
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], label = labels[key], markersize = 12,  
                           mfc = 'None', mew = 2, alpha = 0.5)            
        else:
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], markersize = 12, mew = 0, 
                           label = labels[key], alpha = 0.5)
        axScatter.plot(np.linspace(7.5,12.5), 
                       np.ones(len(np.linspace(7.5,12.5))), 'k-.')
        #axScatter.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
        axScatter.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=15, color='k')
        
        axScatter.text(9.52, -2.2, r'Gas Richness', fontsize=15, color='k')
        axScatter.text(9.52, -2.4, r'Threshold Mass', fontsize=15, color='k')

        axScatter.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
                       np.linspace(10**-2.5,10**1.5), 'k-.')
        axScatter.set_ylim(10**-2.5,10**1.5)
        axScatter.set_xlim(7.5,12.5) 
        axScatter.set_yscale("log")
        axScatter.yaxis.set_major_formatter(ScalarFormatter())
        axScatter.legend(loc='upper right',numpoints = 1, fontsize = 15) # bbox_to_anchor=(1,1),
        #if he2_flag:
            #axScatter.legend(loc=2, bbox_to_anchor=(1,1.15),numpoints = 1)

        axScatter.set_xlabel(r'$\rm \log(M_{stellar}/M_{\odot})$')
        axScatter.set_ylabel(r'$\rm M_{gas}/M_{stellar}$')
axScatter.plot(df.logmstar.loc[jwst], 10**df.logmgas.loc[jwst]/10**df.logmstar.loc[jwst], 
                           '*', color = 'lime',markersize = 12, mew = 0, 
                           label = labels[key])
#for i, txt in enumerate(df.name):
#    axScatter.annotate(txt, (df.logmstar.loc[df.name[i]], 
#                             10**df.logmgas.loc[df.name[i]]/10**df.logmstar.loc[df.name[i]]))
#for i, txt in enumerate(jwst):
#    axScatter.annotate(txt, (df.logmstar.loc[jwst[i]], 
#                             10**df.logmgas.loc[jwst[i]]/10**df.logmstar.loc[jwst[i]]))

if not he2_flag:
    bins = np.arange(7.5,11.5,0.25)
    axHistx.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
             bins= bins, linewidth = 5, label = 'All Galaxies')
    mstars = df.iloc[np.where(flags['defstarform'])[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, alpha = 0.3,
                 linewidth = 5, label = labels['defstarform'], 
                    color = colors['defstarform'])
    mstars = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                    (flags['defliner']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, hatch = '/',
                 linewidth = 5, label = 'Traditional AGN', 
                 color = colors['defagn'])
    mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, hatch = '\\',
                 linewidth = 5, label = 'SF-AGN', 
                 color = colors['sftoagn'])
    axHistx.legend(loc='upper right', fontsize = 15) #bbox_to_anchor=(1,1.15), 
    
axHistx.set_ylabel('Number')
axHistx.set_xlim(axScatter.get_xlim())
axHistx.yaxis.set_major_formatter(ScalarFormatter())
 
#G/S vs M_star w/o histogram
#fig, axScatter = plt.subplots()       
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#        if key == 'defstarform':
#            axScatter.plot(sel.logmstar,
#                           10**sel.logmgas/10**sel.logmstar, 
#                           marker[key], label = labels[key], markersize = 8, 
#                           alpha = 0.3, mew = 0)
#        
#        #elif key == 'sftoagn':
#        #    axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
#        #                   label = labels[key], markersize = 15, mew = 0)
#                
#        elif key == 'agntosf':
#            axScatter.plot(sel.logmstar,
#                           10**sel.logmgas/10**sel.logmstar, 
#                           marker[key], label = labels[key], markersize = 10, 
#                           mew = 1, mec = 'y')            
#        elif key == 'heiisel':
#            axScatter.plot(sel.logmstar,
#                           10**sel.logmgas/10**sel.logmstar, 
#                           marker[key], label = labels[key], markersize = 12,  
#                           mfc = 'None', mew = 2)            
#        else:
#            axScatter.plot(sel.logmstar,
#                           10**sel.logmgas/10**sel.logmstar, 
#                           marker[key], markersize = 12, mew = 0, 
#                           label = labels[key])
#        axScatter.plot(np.linspace(7.5,11.5), 
#                       np.ones(len(np.linspace(7.5,11.5))), 'k-.')
#        #axScatter.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=14, color='k')
#        
#        axScatter.text(9.52, -2.2, r'Gas Richness', fontsize=14, color='k')
#        axScatter.text(9.52, -2.4, r'Threshold Mass', fontsize=14, color='k')
#
#        axScatter.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
#                       np.linspace(10**-2.5,10**1.5), 'k-.')
#        axScatter.set_ylim(10**-2.5,10**1.5)
#        axScatter.set_xlim(7.5,12.5) 
#        axScatter.set_yscale("log")
#        axScatter.yaxis.set_major_formatter(ScalarFormatter())
#        #if he2_flag:
#            #axScatter.legend(loc=2, bbox_to_anchor=(1,1.15),numpoints = 1)
#
#        axScatter.set_xlabel(r'$\rm \log(M_{stellar}/M_\odot)$', fontsize = 22)
#        axScatter.set_ylabel(r'$\rm M_{gas}/M_{stellar}$',fontsize = 22)



#################

lowsfagn = df[flags.sftoagn & (10**df.logmgas/10**df.logmstar < 1.0)]

#FSMGR vs SSFR
#plt.figure()
#
#for key in keys:
#    sel = df.iloc[np.where(flags[key])[0]]
#    ssfr = 10**(np.log10(sel.sfr_nuv) - sel.logmstar)*10**9
#    if key == 'defstarform':
#        plt.plot(ssfr, sel.meanfsmgr,
#             marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
#             color = colors[key], label = labels[key])
#    elif key == 'agntosf': 
#        plt.plot(ssfr, sel.meanfsmgr,
#             marker[key], markersize = 10, mew = 1, color = colors[key],
#             mec = 'y', label = labels[key])
#    else:
#        plt.plot(ssfr, sel.meanfsmgr,
#             marker[key], markersize = 10, mew = 0, color = colors[key],
#             label = labels[key])
#    #plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
#    plt.plot(np.linspace(-12,-8), 1+0*np.linspace(-12,-8), 'k-.')
#    #plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
#    #         rotation = 'vertical')
#    plt.text(-11.9, 1.5, r'Stellar Mass Doubled in last Gyr', 
#             fontsize=14, color='k')
#    
#    plt.xlabel('SSFR (Short Term SF)', size = 22)
#    plt.ylabel('FSMGR (Long Term SF)', size = 22)
#    plt.yscale('log')
#    plt.xscale('log')
#    yticks = plt.yticks()[0]
#    plt.yticks(yticks, np.around(yticks,2))
#    plt.ylim(10**-3, 10**2)
#    plt.xlim(10**-3,2.7)
#newyaxis = np.arange(np.min(df.meanfsmgr)-0.05, np.max(df.meanfsmgr)+0.05,0.01)
#newxaxis = np.ones(len(newyaxis))
#plt.plot(10**-0.26*newxaxis, newyaxis, 'k--', lw = 3)
#    #plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)

ssfr_st = 10**(np.log10(df.sfr_nuv) - df.logmstar)
#ssfr_st = 10**(np.log10(df.sfr_nuv_wise) - df.logmstar)
ssfr_lt = 1/(1+(1/(df.meanfsmgr)))
fsmgr_st_minus1 = 100*(10**6)*(10**-10)/(0.1*1e9*(1-(100*(10**6)*10**-10))) #using ssfr instead of sfr
fsmgr_st = 100*(1e6)*(df.sfr_nuv)/((10**df.logmstar - (100*1e6*df.sfr_nuv)) *0.1*1e9)
fsmgr_lt = df.meanfsmgr
rescat = readsav("resolvecatalog_021622.dat",python_dict=True)
fullfsmgr_st = 100*(1e6)*(rescat['sfr_nuv'])/((rescat['mstars'] - (100*1e6*rescat['sfr_nuv'])) *0.1*1e9)
fullfsmgr_lt = rescat['meanfsmgr']


oldsfr = df['sfr_nuv'].copy()
#Substitute RESOLVE SFR_NUV_WISE to SEL-DF SFR_NUV from database
rescat = readsav("resolvecatalog_031622.dat",python_dict=True)
rescatphot = readsav("resolvecatalogphot_021622.dat",python_dict=True)
match, rescatndx, resselndx = np.intersect1d(rescat['name'], df.name, return_indices=True)
rescatsel_sfr = rescat['sfr_nuv_wise'][rescatndx]
df.loc[df.name.iloc[resselndx],'sfr_nuv'] = rescatsel_sfr
df.loc[df.name.iloc[resselndx],'nuvmag'] = rescatphot['deextrestnuvmag'][rescatndx]


ecocat = pd.read_csv("sfr_nuv_wisew_ecoSEL_newsnr.txt")
ecocat['sfr_nuv_wisew'][(ecocat['mw4w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw4w'] < 0)]
ecocat['sfr_nuv_wisew'][(ecocat['mw3w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw3w'] < 0)]
match, ecocatndx, ecoselndx = np.intersect1d(ecocat['econame'], df.name, return_indices=True)
ecocatsel_sfr = np.array(ecocat['sfr_nuv_wisew'][ecocatndx])
ecocatsel_sfr_err = np.array(ecocat['sfr_nuv_wisew_err'][ecocatndx])
df.loc[df.name.iloc[ecoselndx],'sfr_nuv'] = ecocatsel_sfr

newsfr = df['sfr_nuv'].copy()

ssfr_st = 10**(np.log10(df.sfr_nuv) - df.logmstar)
#ssfr_st = 10**(np.log10(df.sfr_nuv_wise) - df.logmstar)
ssfr_lt = 1/(1+(1/(df.meanfsmgr)))
fsmgr_st_minushalf = 100*(10**6)*(10**-10)/(0.1*1e9*(1-(100*(10**6)*0.5*10**-10))) #using ssfr instead of sfr
fsmgr_st = 100*(1e6)*(df.sfr_nuv)/((10**df.logmstar - (100*1e6*df.sfr_nuv)) *0.1*1e9)
fsmgr_lt = df.meanfsmgr

plt.figure()
full = pd.read_csv('ECO+RESOLVE_inobssample.csv')
fullfsmgr_st = 100*(1e6)*(full['sfr_nuv_wise'])/((10**full['logmstar'] - 
                   (100*1e6*full['sfr_nuv_wise'])) *0.1*1e9)
fullfsmgr_lt = full['meanfsmgr']

diff = (newsfr/oldsfr) > 4
#plt.suptitle(sdsscat+' catalog')
#fsmgr_st = ssfr_st*10**9
#fsmgr_lt = ssfr_lt
#fsmgr_st_minus1 = 10**-0.5
#plt.plot(fullfsmgr_st, fullfsmgr_lt,'.', color = 'g', ms = 10, alpha = 0.5)

for key in keys:
    sel = np.array(df.name.loc[flags[key]])
    plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
            marker[key], color = colors[key], ms = 10, 
            label = labels[key], alpha = 0.7, mew = 0)
#plt.plot(fsmgr_st.loc[diff], fsmgr_lt.loc[diff],
#         'o', color = 'k', ms = 12, mew = 3, mfc = 'none')
#names = ['ECO09516', 'ECO09005', 'ECO02961', 'ECO12402', 'ECO12194']
#names = ['rs0808', 'ECO07935', 'ECO08187', 'ECO12402', 'ECO11902', 'ECO03435']
#names = list(df.name[diff & (fsmgr_lt < 0.005)])
names = ['rs0808', 'ECO07935', 'ECO12402']#, 'ECO11902']
names = ['ECO12402']#, 'ECO11902']
names = ['ECO03435']#, 'ECO11902']
#names = midiragnnames
#for name in names:
#    plt.plot(fsmgr_st.loc[name], fsmgr_lt.loc[name],
#         'o', color = 'black', ms = 12, mew = 3, mfc = 'none')
#
#plt.plot(fsmgr_st.loc[names], fsmgr_lt.loc[names],
#         'o', color = 'green', ms = 12, mew = 3, mfc = 'none')

xaxis = np.arange(np.min(fsmgr_st)-0.05, np.max(fsmgr_st)+0.05,0.01)
yaxis = np.ones(len(xaxis))
plt.plot(xaxis, yaxis, 'k--', lw = 3)
newyaxis = np.arange(np.min(fsmgr_lt)-0.05, np.max(fsmgr_lt)+0.05,0.01)
newxaxis = np.ones(len(newyaxis))
#plt.plot(fsmgr_st_minushalf*newxaxis, newyaxis, 'k--', lw = 3)
plt.xlim(np.min(fsmgr_st)*0.75, np.max(fsmgr_st)*1.1)
plt.ylim(np.min(fsmgr_lt)*0.75, np.max(fsmgr_lt)*1.1)
plt.text(1.5e-12, 1.25, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
plt.tick_params(which = 'minor', length=4, width=2)

#plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
#            '*',color = 'lime', ms = 10, label = labels[key])
plt.xlabel(r'Short Term FSMGR $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr) \times 0.1 Gyr}\right)$')
plt.ylabel(r'Long Term FSMGR $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr) \times 1 Gyr}\right)$')
dwarf = df.name[df.logmstar < 9.5]
#for i, txt in enumerate(df.name):
#    plt.annotate(txt, (fsmgr_st.loc[df.name[i]], 
#                             fsmgr_lt.loc[df.name[i]]))
#for i, txt in enumerate(dwarf):
#    plt.annotate(txt, (fsmgr_st.loc[dwarf[i]], 
#                             fsmgr_lt.loc[dwarf[i]]))

from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mpcolors

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 
plt.figure()
#plt.suptitle('AGNFRAC = 0 only')

ax1 = plt.subplot()
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_spectra/ECO+RESOLVE_snr5_jhu_csf_M5grid_final_LOGZ.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
#Z = pd.read_csv("C:/Users/mugdhapolimera/github/nebulabayes/RESOLVE_jhu_snr5_coincident_open_M-7_0/RESOLVE_jhu_snr5_coincident_open_M-7_0_LOGZ.txt", 
#                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
Z.index = Z.name
Z_pp04 = Z.loc[df.name]['Estimate']+8.76-0.11
#Z_pp04[bad] = -99
ymin, ymax = (7.8, 9.2)
xmin,xmax = (7.5,11.5)

X,Y,Z = density_estimation(df.logmstar.loc[flags.defstarform],
                           Z_pp04.loc[flags.defstarform])
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
for key in keys:
    if key != 'defstarform':
        if key == 'agntosf':
            ax1.plot(df.logmstar.loc[flags[key]], Z_pp04.loc[flags[key]], \
                     marker[key], markersize = 10, label = labels[key], 
                     color = colors[key], alpha = 1.0, mec = 'none')
        else:
            ax1.plot(df.logmstar.loc[flags[key]], Z_pp04.loc[flags[key]], \
                     marker[key], markersize = 10, label = labels[key], 
                     color = colors[key], alpha = 0.5, mec = 'none')

    #ax1.plot(df.logmstar.loc[jwst], Z_pp04.loc[jwst],
#         '*', color = 'lime', markersize = 10, label = labels[key])
ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 4)
yaxis = np.linspace(np.min(Z_pp04)-0.5,np.max(Z_pp04)+0.5)
xaxis = np.linspace(7, 11.5)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
Z04 = np.log10(0.4)+8.76-0.11
plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
plt.plot(xaxis, (8.76-0.11)*np.ones(len(xaxis)), 'k--', linewidth = 2)
#plt.ylim(min(yaxis),max(yaxis))
ax1.set_ylim(7.8,9.2)
plt.xlim(7.5,11.5)
plt.xlabel(r'log(M${\rm_{stellar}}$/M${_\odot}$)')#, fontsize = 22)
plt.ylabel(r'Gas-phase Z (12 + log(O/H))')#, fontsize = 22)
plt.legend(loc='lower right', fontsize = 17)
ax2 = ax1.twinx()
#ax2.set_xticks([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0])
yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
ax2.set_yticks(np.linspace(0,1,len(yticks)))
#ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
float_formatter = lambda x: "%.2f" % x
#xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
N2 = 1.754*yticks - 15.614
N2_label = ["%.2f" % z for z in N2]
ax2.set_yticklabels(N2_label)
ax2.set_ylabel(r'[NII]/H$\alpha$')#, fontsize = 22)
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
#ax1.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
#ax1.plot(m+10,z)

#plt.figure()
#ax1 = plt.subplot()
#N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
#bad = ((N2<-2.5) | (N2>-0.3))
#Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
#Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_dext_jhu_csf_intrinsic_closed_noncoincident_LOGZ.txt", 
#                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
#Z.index = Z.name
#agn = Z.loc[df.name]['Estimate']#+8.76
##Z_pp04[bad] = -99
#for key in keys:
#    if key == 'defstarform':
#        ax1.plot(Z_pp04.loc[flags[key]],agn.loc[flags[key]],
#                 'k.', alpha = 0.3, label = labels[key])
#    else:
#        ax1.plot(Z_pp04.loc[flags[key]],agn.loc[flags[key]],
#         marker[key], markersize = 10, label = labels[key])
#
##yaxis = np.linspace(np.min(Z_pp04)-0.1,np.max(Z_pp04)+0.1)
#xaxis = np.linspace(7.5, 11.5)
##plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
#Z04 = np.log10(0.4)+8.76
##plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
##plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
##plt.ylim(min(yaxis),max(yaxis))
##plt.ylim(7.8,9.2)
##plt.xlim(7.5,11.5)
#plt.ylabel('AGN Fraction', fontsize = 22)
#plt.xlabel(r'Z (12 + log(O/H))', fontsize = 22)
##plt.legend(fontsize = 12)
##ax2 = ax1.twinx()
##ax2.set_xticks([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0])
##yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
##ax2.set_yticks(np.linspace(0,1,len(yticks)))
##ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
##float_formatter = lambda x: "%.2f" % x
##xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
##N2 = 1.754*yticks - 15.614
##N2_label = ["%.2f" % z for z in N2]
##ax2.set_yticklabels(N2_label)
#ax2.set_ylabel(r'[NII]/H$\alpha$', fontsize = 22)

#sfingagn = list(flags.galname[flags.sftoagn])
#Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter_jhu_csf_emergent_open_agn_LOGZ.txt", 
#                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
#Z.index = Z.name
#Z_pp04 = Z.loc[df.name]['Estimate']+8.76
#Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter_jhu_csf_emergent_open_LOGZ.txt", 
#                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
#Z.index = Z.name
#Z_pp04_0 = Z.loc[df.name]['Estimate']+8.76
#zdf = pd.DataFrame(data = {'name' : sfingagn, 
#                           'Z w/ AGNFRAC = 0' : Z_pp04_0.loc[sfingagn], 
#                           'AGNFRAC' : agn.loc[sfingagn],
#                           'Z w/ all AGNFRAC' : Z_pp04.loc[sfingagn]})
#zdf['Z w/ all AGNFRAC (solar units)'] = 10**(zdf['Z w/ all AGNFRAC'] - 8.76)
#zdf['Z w/ AGNFRAC = 0(solar units)'] = 10**(zdf['Z w/ AGNFRAC = 0'] - 8.76)
#
#zdf.to_csv("RESOLVE_SFing-AGN_Metallicities.csv")

#plt.figure()
#n, bins, patches = plt.hist(Z_pp04-8.76, color = 'black', bins = 'fd', 
#                            histtype = 'step', linewidth = 3)
#plt.hist(Z_pp04.loc[flags['sftoagn']]-8.76, color = 'orange', bins = bins,
#         histtype = 'step',hatch = '\\', linewidth = 3)
#yaxis = np.arange(200)
#xaxis = np.ones(len(yaxis))
#plt.plot(np.median(Z_pp04.loc[flags['sftoagn']]-8.76)*xaxis,yaxis, 'k')
#plt.plot(np.median(Z_pp04-8.76)*xaxis,yaxis,'orange')
#plt.yscale('log')

plt.figure()
#plt.suptitle(sdsscat+' catalog')
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
        if key == 'defstarform':        
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 
                     marker[key], markersize = 5, alpha = 0.3, mec = 'none') #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 
                     marker[key], markersize = 10, label = labels[key], alpha = 0.3, mec = 'none') #central
        elif key == 'heiisel':
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], 
                     marker[key], markersize = 5, mfc ='none', mew = 2, alpha = 0.7) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], 
                     marker[key], markersize = 10, mfc ='none', mew = 2,label = labels[key], alpha = 0.5) #central
            
        else:
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], 
                     markersize = 5, alpha = 0.7, mec = 'none') #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], 
                     markersize = 10, label = labels[key], alpha = 0.7, mec = 'none') #central

#        plt.plot(sel.logmh[sel.fc ==0], mbary[sel.fc==0], marker[key], markersize = 5) #other group galaxies
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


#plt.figure()
#ax1 = plt.subplot()
#agnfull = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_spectra/ECO+RESOLVE_snr5_jhu_20Myrssp_M5grid_agn_AGNFRAC.txt", 
#                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
#agnfull.index = agnfull.name
#agn_err = agnfull[["err_up","err_down"]]
#agn = agnfull['Estimate']
##ymin, ymax = (7.8, 9.2)
##xmin,xmax = (7.5,11.5)
##X,Y,Z = density_estimation(df.logmstar.loc[flags.defstarform],
##                           agn.loc[flags.defstarform])
#lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
##def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
##  	new_cmap = mpcolors.LinearSegmentedColormap.from_list(
##        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
##        cmap(np.linspace(minval, maxval, n)))
##  	return new_cmap
##sf_colors_map = truncate_colormap(cm.gray,n=11)
#for key in keys:
#    if key != 'defstarform':
#        if key == 'agntosf':
#            ax1.plot(df.logmstar.loc[flags[key]], agn.loc[flags[key]], \
#                     marker[key], markersize = 10, label = labels[key], 
#                     color = colors[key], alpha = 1.0, mec = 'none')
#        else:
#            ax1.errorbar(df.logmstar.loc[flags[key]], agn.loc[flags[key]], \
#                     yerr = [np.array(agn_err.loc[flags[key]].err_down),\
#                             np.array(agn_err.loc[flags[key]].err_up)],
#                    fmt = marker[key], markersize = 10, label = labels[key], 
#                     color = colors[key], alpha = 0.5, mec = 'none')
#
##ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 4)
#ax1.set_ylim(0,1.1)
##plt.xlim(7.5,11.5)
#plt.xlabel(r'log(M${\rm_{stellar}}$/M${_\odot}$)')#, fontsize = 22)
#plt.ylabel(r'Z (12 + log(O/H))')#, fontsize = 22)
#plt.legend(loc='upper left', fontsize = 20)




xmin = 7.5
xmax = 11.5
ymin = 0
ymax = 3
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_inobssample.csv'
#df = pd.read_csv(inputfile)
#df.index = df.name
#print len(df)
#ra=df.radeg
#dec=df.dedeg
#flinsample = df.fl_insample
#grpcz = df.grpcz
#cz = df.cz
#infall = (ra > 22*15.) | (ra < 3*15.)
#inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
#mgas = df.logmgas
#mstars = df.logmstar
#mbary = 10**mgas + 10**mstars
#inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (np.log10(mbary) > 9.2)
##(((flinsample | (np.log10(mbary) > 9.0)) & infall) | \
##        ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#df = df[inobssample]
#u_r = df['modelu_rcorr']
#
#X,Y,Z = density_estimation(df.logmstar,u_r)
#fig,ax = plt.subplots()#(figsize=(8,8))
#plt.imshow(np.rot90(Z), cmap='bone_r',                                                    
#          extent=[xmin, xmax, ymin, ymax], interpolation='gaussian')
#
#sfagnplt, = ax.plot(df.logmstar.loc[list(flags[flags.sftoagn].galname)], 
#        u_r.loc[list(flags[flags.sftoagn].galname)],
#         'bs', mec = 'k', markersize = 10)#, label = 'SFing-AGN'),
#ax.plot(df.logmstar.loc[names], 
#        u_r.loc[names],
#         'go', mec = 'k', markersize = 10)
#ax.axhline(y = 1.4)
#ax.legend(handles=[sfagnplt], 
#          labels = ['SFing-AGN', 'Targets for this proposal','SFing-AGN in SAMI'],
#          fontsize = 12, loc = 'upper left')
#plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
#plt.ylabel('u-r', fontsize = 15)
#plt.clim(0,1.8)
#plt.contour(X, Y, Z, cmap='summer')
#
#
#ecodat = readsav('eco_wresa_032918.dat')
#match, ecodatndx, dfndx = np.intersect1d(ecodat['econames'], df.name, return_indices = True)
#df.loc[df.name.iloc[dfndx],'nuvmag'] = ecodat['rpdextnuvmagnew'][ecodatndx]
#dfu_r = df.nuvmag - df.absrmag
#
##u_r = ecodat['rpdextnuvmagnew'].ravel() - ecodat['rpdextrmagnew']
#u_r = rescatphot['deextrestnuvmag'] - rescatphot['deextrestrmag']
##u_r = df['nuvmag'] - df['absrmag']
##u_r = u_r[df.nuvmag != 0]
##mstar= ecodat['rpgoodmstarsnew']#eco.logmstar[df.nuvmag != 0]
##mstar= df.logmstar
#mstar= np.log10(rescat['mstars'])
#
#X,Y,Z = density_estimation(mstar,u_r)
#fig,ax = plt.subplots()#(figsize=(8,8))
#plt.imshow(np.rot90(Z), cmap='bone_r',                                                    
#          extent=[xmin, xmax, ymin, ymax], interpolation='gaussian')
#
#sfagnplt, = ax.plot(df.logmstar.loc[list(flags[flags.sftoagn].galname)], 
#        dfu_r.loc[list(flags[flags.sftoagn].galname)],
#         'bs', mec = 'k', markersize = 10)#, label = 'SFing-AGN'),
#ax.plot(df.logmstar.loc[names], 
#        dfu_r.loc[names],
#         'go', mec = 'k', markersize = 10)
#ax.axhline(y = 1.4)
#plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
#plt.ylabel('u-r', fontsize = 15)
#plt.clim(0,1.8)
#plt.contour(X, Y, Z, cmap='summer')
#
#ecowise = readsav('ecoSEL_wise_wisemask_031022.dat')
#ecowise = ecowise.resolve_wise.final[0]
##match, ecowisendx, dfndx = np.intersect1d(ecowise['econame'], df.name, return_indices = True)

ecoir = pd.read_csv('ecoSEL_irphot.txt')
ecoir.index = ecocat.econame

threshold = 5
snr = (ecoir['mw1w']/ecoir['emw1w'] > threshold) & (ecoir['mw2w']/ecoir['emw2w'] > threshold) & \
        (ecoir['mw3w']/ecoir['emw3w'] > threshold) & (ecoir['mw4w']/ecoir['emw4w'] > threshold)
 

w12 = ecoir['mw1w'] - ecoir['mw2w']
w23 = ecoir['mw2w'] - ecoir['mw3w']
w12_err = np.sqrt(ecoir['emw1w']**2 + ecoir['emw2w']**2)
w23_err = np.sqrt(ecoir['emw2w']**2 + ecoir['emw3w']**2)

#mid-IR AGN if data+/-error satisfies the AGN criteria
#Stern OR Jarrett
#midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & (w12 < 1.7) & 
#                            (0.1*w23 + 0.38 < w12)) | 
#           (w12 >= 0.52) & (w12 >= (5.78*w23) -24.50))

#Stern OR Jarrett OR Satyapal
#    midiragn = ((w12-w12_err >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
#                 (w12-w12_err < 1.7) & (0.1*w23 + 0.38 < w12-w12_err)) | \
#               (w12-w12_err >= 0.52) & (w12-w12_err >= (5.78*w23) -24.50))
#Without error consideration
midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
             (w12< 1.7) & (0.1*w23 + 0.38 < w12)) | \
           (w12>= 0.52) & (w12>= (5.78*w23) -24.50))

def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]

def satyapalx(x):
    return 0.52*np.ones(len(x))

def satyapaly(x):
    return 5.78*x -24.50
ax = plt.subplot(111)
ax.plot(xaxis, stern(xaxis), 'k-.')#, label = 'Stern12')
xaxis = np.linspace(min(w23)-0.1,max(w23)+0.1)
#yaxis = np.linspace(min(w12), max(w12))
yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
ax.text(5.75,1.0,'St12', fontsize = 15)
ax.plot(xaxis, satyapalx(xaxis), 'k')#, label = 'Satyapal18')
ax.text(5.75,0.3,'Sa14', fontsize = 15)
ax.text(4.7,2.0 ,'Sa18', fontsize = 15)
xaxis = np.linspace(4.3287,max(w23))
ax.plot(xaxis, satyapaly(xaxis), 'k')

xaxis = np.linspace(2.2,4.2)
ax.plot(jarretx(yaxis)[0], yaxis, 'k--', jarretx(yaxis)[1], yaxis, 'k--')
ax.plot(xaxis, jarrety(xaxis)[0], 'k--')
ax.plot(xaxis, jarrety(xaxis)[1],'k--')#, label = 'Jarrett15')
ax.text(3.5,1.85,'J11', fontsize = 15)
ax.set_xlabel('W2 - W3')
ax.set_ylabel('W1 - W2')
ax.set_ylim(min(w12)-0.1, max(w12)+0.1)
#plt.errorbar(w23,w12,fmt = 'bo', xerr = w23_err,
#             yerr = w12_err, label = 'Galaxies with reliale WISE mags')
#plt.errorbar(w23[midiragn],w12[midiragn],fmt = 'rs', xerr = w23_err[midiragn],
#             yerr = w12_err[midiragn], label = 'Mid-IR AGN')
#plt.errorbar(w23['rs0107'],w12['rs0107'],fmt = 'ks', xerr = w23_err['rs0107'],
#             yerr = w12_err['rs0107'], label = 'rs0107')
#plt.plot((fulldf[~good]['mw2'] - fulldf[~good]['mw3']),
#         (fulldf[~good]['mw1'] - fulldf[~good]['mw2']),
#         'k.', alpha = 0.3, mec = 'none',
#         label = 'Galaxies without reliable WISE mags')
ax.plot(w23,w12,'gp', alpha = 0.3, ms = 10, mec = 'none',
         label = 'Mid-IR SF')
ax.plot(w23[midiragn],w12[midiragn],'p', color = 'orange', ms = 10, mec = 'none',
         label = 'Mid-IR AGN')
dwarfs = df.logmstar < 9.5
dwarfagn = dwarfs & midiragn
ax.plot(w23[dwarfagn],w12[dwarfagn],'kp', ms = 12, mec = 'k', mfc = 'none',
             label = 'Mid-IR Dwarf AGN')

#plt.plot(w23['rs0107'],w12['rs0107'],fmt = 'ks', label = 'rs0107')
ax.set_ylim(-1.2, 2.2)
ax.set_xlim(-1.25,6.5)#min(w23)-0.1,max(w23))

ax.legend(loc = 'lower right', fontsize = 18)


ax.plot(w23[names],w12[names],'o', color = 'black', ms = 10, mec = 'none',
         label = 'Mid-IR AGN')




