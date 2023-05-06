# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:25:30 2018

@author: mugdhapolimera

This code explores the properties of galaxies categorized using BPT plots.
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

#path = os.path.dirname(os.path.realpath(__file__))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
he2_flag = 0
if he2_flag:
    flags = pd.read_csv('resolve_emlineclass_filtered_he2.csv')
else:
    flags = pd.read_csv('resolve_emlineclass_full_snr5_master_new.csv')

flags.index = flags.galname
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_master_new.csv'
full_df = pd.read_csv(inputfile)
full_df.index = full_df.name
df = full_df.loc[flags.galname]

results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_filtered.txt'
full_Z = np.genfromtxt(results_izi, dtype = None, names= True)
Z = full_Z[list(x for x in range(len(full_Z)) if full_Z['Name'][x] in list(flags.galname))]
        

#os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
#results_nb = pd.read_csv("results_izi_prior/RESOLVE_param_estimates.csv")
results_nb = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_master_csf_nicholls_fstar30_agn.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])

#Z_index = (results_nb['Estimate'] == 'LOGZ')
#full_Z_nb = results_nb[Z_index]
results_nb.index = results_nb.name
Z_nb = results_nb.loc[flags.galname]

results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_prior2.txt'
full_Z1 = np.genfromtxt(results_izi, dtype = None, names= True)
Z1 = full_Z1[list(x for x in range(len(full_Z1)) if full_Z1['Name'][x] in list(flags.galname))]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}


#Mass Distribution of all galaxies and BPT classified galaxies
plt.figure()
bins = np.arange(7.25,11.5,0.25)
plt.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
         bins= bins, linewidth = 5, label = 'All Galaxies')
for key in keys:
        mstars = df.iloc[np.where(flags[key])[0]].logmstar
        plt.hist(mstars, histtype = 'step',
                          bins = bins,linewidth = 5, label = key)
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Stellar Mass')
        plt.ylabel('Number')
plt.figure()
bins = np.arange(7,11,0.25)
plt.hist(df.logmgas, histtype = 'stepfilled', alpha = 0.1,
         bins = bins, label = 'All Galaxies')
for key in keys:
        mgas = df.iloc[np.where(flags[key])[0]].logmgas
        plt.hist(mgas, histtype = 'step', 
                         bins = bins,linewidth = 5, label = key)
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Gas Mass')
        plt.ylabel('Number')

bins = np.arange(0.775,1.2,0.025)
plt.figure()
plt.hist(df.logmgas/df.logmstar, bins = bins, alpha = 0.1, 
         histtype = 'stepfilled', label = 'All Galaxies')
for key in keys:
        mgas = df.iloc[np.where(flags[key])[0]].logmgas
        mstars = df.iloc[np.where(flags[key])[0]].logmstar
        
        plt.hist(mgas/mstars, histtype = 'step', bins = bins,
                         linewidth = 5, label = key,)
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Gas/Stellar Mass Ratio')
        plt.ylabel('Number')

        
plt.figure()
bins = np.arange(8.25,11.5,0.25)
plt.hist(np.log10(10**df.logmstar + 10**df.logmgas), histtype = 'stepfilled', 
         alpha = 0.1, bins = bins, label = 'All Galaxies')
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        mbary = np.log10(10**sel.logmstar + 10** sel.logmgas)
        plt.hist(mbary, histtype = 'step', 
                         bins = bins,linewidth = 5, label = key)
        plt.legend()
        plt.yscale ('log')
        plt.xlabel('Baryonic Mass')
        plt.ylabel('Number')
        
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            plt.plot(sel.logmstar,sel.logmgas, marker[key], label = key, 
                     alpha = 0.3)
        else:
            plt.plot(sel.logmstar,sel.logmgas, marker[key],
                     label = key)
        plt.plot(np.linspace(7.5,11.5), np.linspace(7.5,11.5), 'k-.')
        if he2_flag:
            plt.legend(loc = 4, numpoints = 1)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Gas Mass')

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(9)

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx, yscale = 'log')
axHisty = plt.axes(rect_histy, xscale = 'log')

# no labels
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], label = key, 
                     alpha = 0.3)
        elif key == 'heiisel':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], label = key, 
                     mfc = 'None', mew = 2)            
        else:
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key],
                     label = key)
        axScatter.plot(np.linspace(7.5,11.5), np.linspace(7.5,11.5), 'k-.')
        if he2_flag:
            axScatter.legend(loc = 4, numpoints = 1)
        axScatter.set_xlabel('Stellar Mass')
        axScatter.set_ylabel('Gas Mass')

axHistx.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
         bins= bins, linewidth = 5, label = 'All Galaxies')
for key in keys:
        mstars = df.iloc[np.where(flags[key])[0]].logmstar
        axHistx.hist(mstars, histtype = 'step',
                          bins = bins,linewidth = 5, label = key)
        #axHistx.legend()
        if he2_flag:
            axHistx.legend(loc=2, bbox_to_anchor=(1,1.15))
        #axHistx.set_xlabel('Stellar Mass')
        axHistx.set_ylabel('Number')

bins = np.arange(7.5,11.5,0.25)
axHisty.hist(df.logmgas, histtype = 'stepfilled', alpha = 0.1,
         bins = bins, label = 'All Galaxies', orientation='horizontal')
for key in keys:
        mgas = df.iloc[np.where(flags[key])[0]].logmgas
        axHisty.hist(mgas, histtype = 'step',orientation='horizontal',
                         bins = bins,linewidth = 5, label = key)
        #axHisty.set_xlabel('Gas Mass')
        axHisty.set_xlabel('Number')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

inputfile = 'C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl'
full = pd.read_pickle(inputfile)
ra=full.radeg
dec=full.dedeg
if 'fl_insample' in full.keys():
    flinsample = full.fl_insample
else:
    flinsample = np.ones(len(full), dtype = bool)
grpcz = full.grpcz
cz = full.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = full.logmgas
mstars = full.logmstar
mbary = 10**mgas + 10**mstars

inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
full = full[inobssample]
full = full[~np.isnan(full.h_alpha_flux_ext)]
full = full[~np.isnan(full.oiii_5007_flux_ext)]
full = full[~np.isnan(full.nii_6584_flux_ext)]
full = full[~np.isnan(full.nii_6548_flux_ext)]
full = full[~np.isnan(full.h_beta_flux_ext)]
full = full[~np.isnan(full.oi_6300_flux_ext)]
full = full[~np.isnan(full.sii_6717_flux_ext)]
full = full[~np.isnan(full.sii_6731_flux_ext)]
ndx = [x for x in full.NAME if x not in list(flags.galname)]
full = full.loc[ndx]
plt.figure(10)
mbary = np.log10(10**full.logmstar + 10**full.logmgas)
if not he2_flag:
    plt.plot(full.logmh[full.fc == 0], mbary[full.fc == 0], 'k+', alpha = 0.3, markersize = 5) #other group galaxies
    plt.plot(full.logmh[full.fc == 1], mbary[full.fc == 1], 'k+', alpha = 0.3, markersize = 10, label = 'Not SF or AGN (Dead)') #central

for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
        if key == 'defstarform':        
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5, alpha = 0.3) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, label = key, alpha = 0.3) #central
        elif key == 'heiisel':
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5, mfc ='none', mew = 2) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, mfc ='none', mew = 2,label = key) #central
            
        else:
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, label = key) #central

#        plt.plot(sel.logmh[sel.fc ==0], mbary[sel.fc==0], marker[key], markersize = 5) #other group galaxies
                
        if he2_flag:
            plt.legend(loc = 2, numpoints = 1)
        plt.xlabel('Group Halo Mass')
        plt.ylabel('Galaxy Baryonic Mass')
        
#plt.figure()
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#            
#        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
#        
#        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
#        plt.legend(borderaxespad=0., loc = 4)
#        plt.xlabel('Stellar Mass')
#        plt.ylabel('Metallicity')
#        plt.title('M-Z Relation using NebulaBayes + Levesque Grid')

#m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
#plt.plot(m_z[:,0], m_z[:,1],'r')
##From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
#m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
#z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
##plt.ylim(-0.2, 1.2)
#plt.plot(m+10,z)
#
#plt.figure()
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#            
#        Z_sel_ndx = Z[list(x for x in range(len(Z)) if Z['Name'][x] in sel.NAME)]
#        Z_sel = Z_sel_ndx['Z_Estimate']
#        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
#        plt.legend(borderaxespad=0., loc = 4)
#        plt.xlabel('Stellar Mass')
#        plt.ylabel('Metallicity')
#        plt.title('M-Z Relation using IZI (Python + GPy) without Prior')
#m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
#plt.plot(m_z[:,0], m_z[:,1],'r')
##From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
#m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
#z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(8, 9.2)
#plt.plot(m+10,z)

#plt.figure()
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#            
#        Z_sel_ndx = Z1[list(x for x in range(len(Z1)) if Z1['Name'][x] in sel.NAME)]
#        Z_sel = Z_sel_ndx['Z_Estimate']
#        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
#        plt.legend(borderaxespad=0., loc = 4)
#        plt.xlabel('Stellar Mass')
#        plt.ylabel('Metallicity')
#        plt.title('M-Z Relation using IZI (Python + GPy) with Prior')
#m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
#plt.plot(m_z[:,0], m_z[:,1],'r')
##From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
#m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
#z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(8, 9.2)
#plt.plot(m+10,z)
#
#plt.figure()
#bins = np.arange(7.6,9.1,0.1)
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
#        plt.hist(Z_sel, histtype = 'step', 
#                         bins = bins,linewidth = 5, label = key)
#        plt.legend()
#        plt.xlabel('Metallicity (NebulaBayes + Levesque Grid)')
#        plt.ylabel('Number')
#
#plt.figure()
#bins = np.arange(7.6,9.1,0.1)
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#            
#        Z_sel_ndx = Z[list(x for x in range(len(Z)) if Z['Name'][x] in sel.NAME)]
#        Z_sel = Z_sel_ndx['Z_Estimate']
#        plt.hist(Z_sel, histtype = 'step', 
#                        bins = bins, linewidth = 5, label = key)
#        plt.legend(borderaxespad=0., loc = 2)
#        plt.xlabel('Metallicity (IZI - Python + Gpy)')
#        plt.ylabel('Number')
#
#plt.figure()
#bins = np.arange(7.6,9.1,0.1)
#for key in keys:
#        sel = df.iloc[np.where(flags[key])[0]]
#            
#        Z_sel_ndx = Z1[list(x for x in range(len(Z1)) if Z1['Name'][x] in sel.NAME)]
#        Z_sel = Z_sel_ndx['Z_Estimate']
#        plt.hist(Z_sel, histtype = 'step', 
#                        bins = bins, linewidth = 5, label = key)
#        plt.legend(borderaxespad=0., loc = 2)
#        plt.xlabel('Metallicity (IZI - Python + Gpy) without Prior')
#        plt.ylabel('Number')

'''plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        #mbary = np.log10(10**sel.logmstar + 10** sel.logmgas)
        bins = np.arange(-6,11,1)
        plt.hist(sel.morph[(sel.morph > -999.0) & (sel.morph < 13.0)], histtype = 'step', 
                         bins= bins,linewidth = 5, label = key)
        plt.legend(loc = 2)
        plt.xlabel('Morphology')
        plt.ylabel('Number')

plt.figure()
bins = np.arange(0,1.5,0.5)
print bins
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        morphel = sel.morphel
        morphel[morphel == 'E'] = 0
        morphel[morphel == 'L'] = 1
        plt.hist(morphel, histtype = 'step', 
                       bins = bins, linewidth = 5, label = key)
        plt.legend(loc = 2)
        plt.xlabel('Morphology')
        plt.ylabel('Number')
'''