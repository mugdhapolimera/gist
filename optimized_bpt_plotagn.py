# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:51:04 2022

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 17:31:14 2021

@author: mugdhapolimera
"""

#This program makes a Line-Ratio diagram (also known as a BPT plot or Kewley diagram)
#with labels using data from the RESOLVE survey to classify galaxies as LINERs,
#Seyferts, Composites, or AGNs on the basis of their flux ratio for distinction.

#Original code from Ashley Bittner 08/03/2017
#Edited version from Margie Bruff 01/07/2018
#Updated by Carlynn Ferguson 03/30/2018

#suggested use of python debugger to understand the code more thoroughly
#KDE plotting 
#https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density/
#import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import sys
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})
from astropy.stats import binom_conf_interval

#define demarcation function: log_NII_HA vs. log_OIII_HB
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2halocus(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.1 + (0.61 / (log_NII_HA + 0.08))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
#    return 0.57 + (0.13 / (log_NII_HA - 0.003))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
    #return 0.58 + (0.04 / (log_SII_HA +0.012))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
    #return 0.61 + (0.056 / (log_OI_HA + 0.40))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59

def he2hbmain(log_NII_HA):
    return -1.22+1.0/(8.92*log_NII_HA+1.32)

def he2hbmainclass(log_NII_HA):
    refn2ha = np.linspace(-3.0, -0.15)
    main = he2hbmain(refn2ha)
    line = np.poly1d(np.polyfit(refn2ha, main, 15))
    return line(log_NII_HA) 

def he2hblimit(log_NII_HA):
        return -1.07+1.0/(8.92*log_NII_HA-0.95)

def he2hblimitclass(log_NII_HA):
    refn2ha = np.linspace(-3.0, -0.15)
    limit = he2hblimit(refn2ha)
    line = np.poly1d(np.polyfit(refn2ha, limit, 15))
    return line(log_NII_HA) 

def ratioerror(num,num_err,den, den_err):
    err_num2 = (num_err/(num*np.log(10)))**2
    err_den2 = (den_err/(den*np.log(10)))**2
    return np.sqrt(err_num2 + err_den2)
    # note that the ratio uses only the stronger line, but for S/N reasons we add
    # the weaker and multiply by 3/4 since Chris Richardson says the canonical
    # line ratio is 3:1 (this needs to be updated with a more precise number)
def bpt_plots(inputfile, outputfile,  selinputfile, s06outputfile,  midirfile, 
              eco, resolve, full, sdsscat, save, ax1, ax2, ax3, simple):
    df = pd.read_csv(inputfile)
        #define alternate catalog names
    if 'name' in df.keys():
        df['NAME'] = df['name']
    if 'CATID' in df.keys():
        df['NAME'] = df['CATID']
    if 'NAME' in df.keys():
        df['NAME'] = df['name']

    if 'name' in df.keys():
        name = df['name']
    elif 'jhu_index' in df.keys():
        df['name']= df['jhu_index']
    df.index = df.name

    s06 = pd.read_csv(s06outputfile)
    s06.index = s06.galname
#    s06agn = np.array(s06[~s06.defstarform].galname)
    s06agn = np.array(s06.galname[s06.composite])

    if midirfile != '':
        midir = pd.read_csv(midirfile)
        midir.index = midir.name
        midiragn = np.array(midir.name.iloc[np.where(midir.agnflag == True)[0]])
    else:
        midiragn = []

    if sdsscat == 'port':
        nii = df['Flux_NII_6583']
        nii_sum = (df['Flux_NII_6583']+ df['Flux_NII_6547'])*3./4
        nii_sum_err = (np.sqrt(df['Flux_NII_6547_Err']**2 + df['Flux_NII_6583_Err']**2))*3./4
        oiii = df['Flux_OIII_5006']
        oiii_err = df['Flux_OIII_5006_Err']
        h_alpha = df['Flux_Ha_6562']
        h_alpha_err = df['Flux_Ha_6562_Err']
        h_beta = df['Flux_Hb_4861']
        h_beta_err = df['Flux_Hb_4861_Err']
        oi = df['Flux_OI_6300']
        oi_err = df['Flux_OI_6300_Err']
        sii_sum = df['Flux_SII_6716'] + df['Flux_SII_6730']
        sii_sum_err = np.sqrt(df['Flux_SII_6716_Err']**2 + df['Flux_SII_6730_Err']**2)
    if sdsscat == 'jhu' or sdsscat == 'nsa' or sdsscat == 'master':
        
        nii = df['nii_6584_flux']
        if 'nii_6548_flux' in df.keys():
            nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
            nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
        else:        
            nii_sum = df['nii_6584_flux']
            nii_sum_err = df['nii_6584_flux_err']

        oiii = df['oiii_5007_flux']
        oiii_err = df['oiii_5007_flux_err']
        h_alpha = df['h_alpha_flux']
        h_alpha_err = df['h_alpha_flux_err']
        h_beta = df['h_beta_flux']
        h_beta_err = df['h_beta_flux_err']
        oi = df['oi_6300_flux']
        oi_err = df['oi_6300_flux_err']
        if 'sii_6717_flux' in df.keys():
            sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
    
            sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
        else:
            sii_sum = df['sii_6731_flux']
        
            sii_sum_err = df['sii_6731_flux_err']
    
    
    data = df.name > 0
        
    if eco: 
        resname = df['resname'] #for eco
        resname = resname != 'notinresolve'
    if resolve:
        econame = df['NAME']#df['econame'] #for resolve
        econame = df['NAME']#econame != 'notineco'

    #print total points shared with alternate catalog
    if full: 
        sel = (np.where(data)[0]) #for eco
    elif eco: 
        sel = (np.where(data & resname)[0]) #for eco
    elif resolve:
        sel = (np.where(data & econame)[0]) #for resolve
    else:
        sel = np.where(data)[0]
    print ''
    print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)
    #df = df[data]
    
    nii = nii[data]
    nii_sum = nii_sum[data]
    nii_sum_err = nii_sum_err[data]
    oiii = oiii[data]
    oiii_err = oiii_err[data]
    oi = oi[data]
    oi_err = oi_err[data]
    sii_sum = sii_sum[data]
    sii_sum_err = sii_sum_err[data]
    h_beta = h_beta[data]
    h_beta_err = h_beta_err[data]
    h_alpha = h_alpha[data]
    h_alpha_err = h_alpha_err[data]
    
    if sdsscat =='nsa':
        df = df[data]
    
        
    sel = pd.read_csv(selinputfile)
    if 'jhu_index' in sel.keys():
        sel['name'] = sel['jhu_index']
    sel.index = sel.name
    onlybpt = set(df.name) - set(sel.name)
    selname = set(sel.name)
#    oi.loc[onlybpt] = np.nan
#    oi_err.loc[onlybpt] = np.nan
#    sii_sum.loc[onlybpt] = np.nan
#    sii_sum_err.loc[onlybpt] = np.nan

#    siisel = pd.read_csv('ECO+RESOLVE_bpt1snr5_siisnr_jhu.csv')
#    siisel.index=siisel.name
#    siisel = siisel['snr']
#    oisel = pd.read_csv('ECO+RESOLVE_bpt1snr5_oisnr_jhu.csv')
#    oisel.index=oisel.name
#    oisel = oisel['snr']
    siisel = np.in1d(df.name[data], sel.name)
    oisel = np.in1d(df.name[data], sel.name)
    
    #bpt = pd.read_csv('ECO+RESOLVE_bpt1snr5_dext_'+sdsscat+'.csv')
    #bpt.index = bpt.name
    #onlybpt = set(bpt.name)
    #selname = set(df.name) - set(bpt.name)
    

    
    #length of data to be used for debugging
    datalen = np.sum(data)
    subsetname = df.name[data]
    
    # data ratios
    #n2ha = np.log10(nii_sum/h_alpha)
    o3hb = np.log10(oiii/h_beta) # always the y-axis
    o1ha = np.log10(oi/h_alpha)
    s2ha = np.log10(sii_sum/h_alpha)
    n2ha = np.log10(nii_sum/h_alpha)
    
    #Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
    #and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.
    
    #NII plot selectors
    compsel1 = (o3hb >= n2hacompmin(n2ha)) & (o3hb <= n2hamain(n2ha))
    sfsel1 = (o3hb < n2hacompmin(n2ha)) & (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
    agnsel1= (o3hb > n2hamain(n2ha))
    #plt.hist(o1ha_err[o1ha_err < 1e5], bins = 'fd')

    #SII plot selectors
    sfsel2 = (o3hb <= s2hamain(s2ha)) & ~compsel1 & siisel
    seyfsel2 = ((o3hb > s2hamain(s2ha)) & (o3hb >= s2halinseyf(s2ha))) & siisel
    linersel2 = ((o3hb > s2hamain(s2ha)) & (o3hb < s2halinseyf(s2ha))) & siisel
    agnsel2 = (o3hb > s2hamain(s2ha)) & ~compsel1 & siisel
    
    #OI plot selectors
    sfsel3 = (o3hb <= o1hamain(o1ha)) & (o1ha < -0.7) & ~compsel1 & oisel
    seyfsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb >= o1halinseyf(o1ha)) & oisel
    linersel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb < o1halinseyf(o1ha)) & oisel
    agnsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & ~compsel1 & oisel
    
    #REFERENCE for cumulative plot selectors
    seyfselr = seyfsel2 & seyfsel3
    linerselr = linersel2 & linersel3
    
    #cumulative plot selectors
    sfsel = (sfsel1 & sfsel2 & sfsel3) | (sfsel1 & ~(siisel | oisel)) 
    compsel = compsel1  #composite galaxies
    seyfsel = agnsel1 & seyfselr #Seyfert AGN galaxies
    linersel = agnsel1 & linerselr #LINER AGN galaxies
    ambigsel1 = sfsel1 & (agnsel2 | agnsel3) #SF in first plot, AGN in subsequent plot
    ambigsel2 = np.array(agnsel1) & (np.array(sfsel2) | np.array(sfsel3)) #AGN in first plot, SF in subsequent plot
    ambagnsel = agnsel1 & ~seyfselr & ~linerselr & ~(sfsel2 | sfsel3) #Ambiguous AGN
    
    sftoagn1 = sfsel1 & agnsel2
    sftoagn2 = sfsel1 & agnsel3
    
    #Save the BPT flags to a CSV file
    emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
    defagn = seyfsel | linersel | ambagnsel
    convagn = compsel | defagn
    
    flags = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                              'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                              'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
                              'sftoagn1':sftoagn1, 'sftoagn2': sftoagn2, 'convagn': convagn})
    if save:
        flags.to_csv(outputfile ,index=False)
#    keys = ['defagn', 'composite', 'sftoagn', 'agntosf','dwarfagn']
    keys = ['convagn', 'sftoagn']
    marker = {'convagn': '1', 'agntosf': 'c^', 'ambigagn': 'ms', 'composite': 'ms', 
              'defagn': 'ro', 
              'defliner': 'yo', 'defseyf': 'co', 'dwarfagn': 'ks',
              'sftoagn': 'bs'}
    markersize = {'convagn': 10, 'agntosf': 14, 'ambigagn': 8, 'composite': 8, 'defagn': 8, 
              'defliner': 8, 'defseyf': 8, 'sftoagn': 8, 
              'dwarfagn':16}
    alpha = {'convagn': 1, 'agntosf': 1, 'ambigagn': 0.5, 'composite': 0.5, 'defagn': 0.5, 
              'defliner': 0.5, 'defseyf': 0.5, 'sftoagn': 0.7, 
              'dwarfagn':1}
    
    markercolors = {'convagn': 'r', 'agntosf': 'c', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
              'defliner': 'y', 'defseyf': 'c', 'dwarfagn': 'k', 'sftoagn': 'b'}
    
    labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
              'composite': 'Composite', 'defagn': 'Traditional AGN', 
              'defliner': 'LINER', 'defseyf': 'Seyfert', 
              'dwarfagn': 'SEL Dwarf AGN', 'defstarform': 'Definite SF', 
              'sftoagn': 'SF-AGN', 'convagn': 'Conventional AGN'}
    
    
    #checking that plotted points are within the total data range
    print ''
    sfselpts = (len(np.where(sfsel)[0]))
    seyfselpts = (len(np.where(seyfsel)[0]))
    linerselpts = (len(np.where(linersel)[0]))
    compselpts = (len(np.where(compsel)[0]))
    agnselpts = (len(np.where(ambagnsel)[0]))
    ambigsel1pts = (len(np.where(ambigsel1)[0]))
    ambigsel2pts = (len(np.where(ambigsel2)[0]))
    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+\
    ambigsel1pts+ambigsel2pts


    sfpercent = float(sfselpts)/float(datalen)*100
    seyfpercent = float(seyfselpts)/float(datalen)*100
    linerpercent = float(linerselpts)/float(datalen)*100
    comppercent = float(compselpts)/float(datalen)*100
    agnpercent = float(agnselpts)/float(datalen)*100
    ambig1percent = float(ambigsel1pts)/float(datalen)*100
    ambig2percent = float(ambigsel2pts)/float(datalen)*100
    #df.index = fulldf.name
    if 'logmstar' in df.keys():
        dwarf = (df.logmstar < 9.5)
        giant = (df.logmstar > 9.5)
    
    agn = (ambigsel1|seyfsel|linersel|ambagnsel|compsel|ambigsel2)
    dwarfagn = dwarf & agn
    giantagn = giant & agn

    s06dwarfagn = np.intersect1d(s06agn, df.name[dwarf])
    midirdwarfagn = np.intersect1d(midiragn, df.name[dwarf])
    
    print ("DATA POINTS: "),datalen
    print ("TOTAL PLOTTED POINTS: "), totalselpts
    print ("TOTAL PLOTTED POINTS OMITTED: "), datalen-totalselpts
    print "* IF NUMBERS ABOVE ARE AT ALL NEGATIVE THEN THERE IS OVERPLOTTING"
    print ("Definite Star Forming: "),sfselpts,("("),round(sfpercent, 2),("%"),(")")
    print ("Composite: "),compselpts, ("("),round(comppercent, 2),("%"),(")")
    print ("SF --> AGN: "), ambigsel1pts, ("("),round(ambig1percent, 2),("%"),(")")
    print ("AGN --> SF: "), ambigsel2pts, ("("),round(ambig2percent, 2),("%"),(")")
    print ("Ambiguous AGN: "),agnselpts, ("("),round(agnpercent, 2),("%"),(")")
    print ("Seyfert: "),seyfselpts, ("("),round(seyfpercent, 2),("%"),(")")
    print ("LINER: "),linerselpts, ("("),round(linerpercent, 2),("%"),(")")
    print ("TOTAL KNOWN AGN: "),linerselpts+seyfselpts+agnselpts, ("("), \
    round(linerpercent+seyfpercent+agnpercent, 2), ("% )")
    print ("POSSIBLE TOTAL AGN: "),linerselpts+seyfselpts+agnselpts+ambigsel1pts+ambigsel2pts,("("),\
    round(linerpercent+seyfpercent+agnpercent+ambig1percent+ambig2percent, 2), ("% )")
    print ("Percent Omitted: "), round((100-(sfpercent+seyfpercent+linerpercent+comppercent+agnpercent+ambig1percent+ambig2percent)), 2), ("%")
    print ''
    
    print ("AGN in Dwarf Galaxies: "), 100*round(np.sum(dwarfagn)/float(np.sum(dwarf)),2), ("%")
    print ("AGN in Giant Galaxies: "), 100*round(np.sum(giantagn)/float(np.sum(giant)),2), ("%")
    print ("AGN in dwarfs: "), np.sum(agn & dwarf)
    print ("Number of Dwarfs:"), np.sum(dwarf)
    
    ###PLOTS###
    #reference points in x-direction for demarcation lines on plots
    if simple:
        refn2ha = np.linspace(-3.0, 0.35)
        refoiha = np.linspace(-2.5, -0.4)
        refsiiha = np.linspace(-2, 0.3,100)
        
        def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
          	new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
          	return new_cmap
        sf_colors_map = truncate_colormap(cm.gray_r)
        ndx = []#'ECO03494']#np.where((df.NAME == 'rs1105') | (df.NAME == 'rs1375'))[0]
        xmin = refn2ha.min(); xmax = refn2ha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        
        #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
        
        #NII/OIII plot
        definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
        ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)

        ax1.set_xlim(-1.5,0.5)
        ax1.set_ylim(-1.0,1.0)
        ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
        ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
        
        for key in keys:
            n2ha_sel = n2ha.loc[flags[key]]
            o3hb_sel = o3hb.loc[flags[key]]
            ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
                     alpha = alpha[key], markersize = markersize[key], 
                     mew = 2, mec = markercolors[key], label = labels[key])

        n2ha_sel = n2ha[dwarfagn]
        o3hb_sel = o3hb[dwarfagn]
        ax1.plot(n2ha_sel, o3hb_sel, 'ks', mfc = 'none',
                 markersize = 14, mew = 2, label = 'Dwarf AGN')

        n2ha_sel = n2ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax1.plot(n2ha_sel, o3hb_sel, 'v',color = 'lime',  mfc = 'none',
                 markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')

        n2ha_sel = n2ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax1.plot(n2ha_sel, o3hb_sel, 'p', color = 'orange', zorder = 10, mfc = 'none',
                 markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
        
        ax1.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
        
        main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Theoretical Maximum Starburst Line (Ke01)')
        composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
        'k--', label = 'Composite Line (Ka03)')
        
        #SII plot
        xmin = refsiiha.min(); xmax = refsiiha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        s2ha = s2ha.loc[siisel]
        o3hb_orig = o3hb.copy()
        o3hb = o3hb.loc[siisel]
        
        s2ha = s2ha[~np.isnan(s2ha)]
        o3hb = o3hb[~np.isnan(s2ha)]
        
        definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
        agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                        ymin_agn:ymax_agn:nbins*1j]
        k = kde.gaussian_kde(agn_contour.T)
        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
        ax2.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map)
        main1, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
        ax2.set_xlim(-1.5,0.5)
        ax2.set_ylim(-1,1)
        ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
        ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                          'k-.', label = 'Liner/Seyfert Division')
        for key in keys:
            s2ha_sel = s2ha.loc[flags[key]]
            o3hb_sel = o3hb.loc[flags[key]]
            ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
                     alpha = alpha[key], markersize = markersize[key], 
                     mew = 2, mec = markercolors[key], label = labels[key])

        s2ha_sel = s2ha[dwarfagn]
        o3hb_sel = o3hb[dwarfagn]
        ax2.plot(s2ha_sel, o3hb_sel, 'ks', mfc = 'none',
                 markersize = 14, mew = 2, label = 'Dwarf AGN')

        s2ha_sel = s2ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax2.plot(s2ha_sel, o3hb_sel, 'v', color = 'lime', mfc = 'none',
                 markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')

        s2ha_sel = s2ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax2.plot(s2ha_sel, o3hb_sel, 'p', mfc = 'none', color = 'orange', zorder = 10, 
                 markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
        
        #OI Plot
        xmin = refoiha.min(); xmax = 0#refoiha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        
        o1ha = o1ha.loc[oisel]
        o3hb = o3hb_orig.loc[oisel]

        o1ha = o1ha[~np.isnan(o1ha) & np.isfinite(o1ha)]
        o3hb = o3hb[~np.isnan(o1ha) & np.isfinite(o1ha)]
        
        definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))

        ax3.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map)
        main1, = ax3.plot(refoiha[refoiha < -0.75], o1hamain(refoiha[refoiha < -0.75]), 'k', 
                          label = 'Ke01 Maximum Starburst Line')
        ax3.set_xlim(-2.0,-0.4)
        ax3.set_ylim(-1,1)
        ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                                       'k-.', label = 'Ke06 Liner/Seyfert Division Line')
        ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)

        for key in keys:
            o1ha_sel = o1ha.loc[flags[key]]
            o3hb_sel = o3hb.loc[flags[key]]
            ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
                     alpha = alpha[key], markersize = markersize[key], 
                     mew = 2, mec = markercolors[key], label = labels[key])

        o1ha_sel = o1ha[dwarfagn]
        o3hb_sel = o3hb[dwarfagn]
        ax3.plot(o1ha_sel, o3hb_sel, 'ks', mfc = 'none',
                 markersize = 14, mew = 2, label = 'Dwarf AGN')

        o1ha_sel = o1ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax3.plot(o1ha_sel, o3hb_sel, 'v', color = 'lime', mfc = 'none',
                 markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')

        o1ha_sel = o1ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax3.plot(o1ha_sel, o3hb_sel, 'p', mfc = 'none', color = 'orange', zorder = 10, 
                 markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
        
        print(100.0*np.sum(dwarfagn)/np.sum(dwarf), \
              100.0*binom_conf_interval(np.sum(dwarfagn),np.sum(dwarf)) - \
              (100.0*np.sum(dwarfagn)/np.sum(dwarf)))
        return (ax1, ax2, ax3)
        
    else:
        refn2ha = np.linspace(-3.0, 0.35)
        refoiha = np.linspace(-2.5, -0.4)
        refsiiha = np.linspace(-2, 0.3,100)
        
        def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
          	new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
          	return new_cmap
        sf_colors_map = truncate_colormap(cm.gray_r)
        ndx = []#'ECO03494']#np.where((df.NAME == 'rs1105') | (df.NAME == 'rs1375'))[0]
        xmin = refn2ha.min(); xmax = refn2ha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        
        #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
        
        #NII/OIII plot
        definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))

        ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
        ax1.set_xlim(-1.5,0.5)
        ax1.set_ylim(-1.0,1.0)
        ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
        ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
        
        for key in keys:
            if key != 'dwarfagn':
                n2ha_sel = n2ha.loc[flags[key]]
                o3hb_sel = o3hb.loc[flags[key]]
                ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
                         markersize = markersize[key], mew = 0, label = labels[key])
            else:
                n2ha_sel = n2ha[dwarfagn]
                o3hb_sel = o3hb[dwarfagn]
                ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
                         markersize = markersize[key], mew = 2, label = labels[key])
        n2ha_sel = n2ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax1.plot(n2ha_sel, o3hb_sel, 'kv', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
        n2ha_sel = n2ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax1.plot(n2ha_sel, o3hb_sel, 'kp', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
        
        ax1.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
        
        main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Theoretical Maximum Starburst Line (Ke01)')
        composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
        'k--', label = 'Composite Line (Ka03)')
        
        #SII plot
        xmin = refsiiha.min(); xmax = refsiiha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        
        definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
        agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                        ymin_agn:ymax_agn:nbins*1j]
        k = kde.gaussian_kde(agn_contour.T)
        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
        ax2.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map)
        main1, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
        ax2.set_xlim(-1.5,0.5)
        ax2.set_ylim(-1,1)
        ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
        ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                          'k-.', label = 'Liner/Seyfert Division')
        for key in keys:
            if key != 'dwarfagn':
                s2ha_sel = s2ha.loc[flags[key]]
                o3hb_sel = o3hb.loc[flags[key]]
                ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
                         markersize = markersize[key], mew = 0, label = labels[key])
            else:
                s2ha_sel = s2ha[dwarfagn]
                o3hb_sel = o3hb[dwarfagn]
                ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
                         markersize = markersize[key], mew = 2, label = labels[key])
        s2ha_sel = s2ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax2.plot(s2ha_sel, o3hb_sel, 'kv', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
        s2ha_sel = s2ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax2.plot(s2ha_sel, o3hb_sel, 'kp', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
        
        #OI Plot
        xmin = refoiha.min(); xmax = 0#refoiha.max()
        ymin = -1.25; ymax = 1.5
        nbins = 50
        
        definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
        k2 = kde.gaussian_kde(definite.T)
        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
        agn_contour = np.column_stack((o1ha[defagn], o3hb[defagn]))
        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                        ymin_agn:ymax_agn:nbins*1j]
        k = kde.gaussian_kde(agn_contour.T)
        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
        ax3.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=sf_colors_map)
        main1, = ax3.plot(refoiha[refoiha < -0.75], o1hamain(refoiha[refoiha < -0.75]), 'k', 
                          label = 'Ke01 Maximum Starburst Line')
        ax3.set_xlim(-2.0,-0.4)
        ax3.set_ylim(-1,1)
        ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                                       'k-.', label = 'Ke06 Liner/Seyfert Division Line')
        ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
        for key in keys:
            if key != 'dwarfagn':
                o1ha_sel = o1ha.loc[flags[key]]
                o3hb_sel = o3hb.loc[flags[key]]
                ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
                         markersize = markersize[key], mew = 0, label = labels[key])
            else:
                o1ha_sel = o1ha[dwarfagn]
                o3hb_sel = o3hb[dwarfagn]
                ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
                         markersize = markersize[key], mew = 2, label = labels[key])
        o1ha_sel = o1ha[s06dwarfagn]
        o3hb_sel = o3hb[s06dwarfagn]
        ax3.plot(o1ha_sel, o3hb_sel, 'kv', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
        o1ha_sel = o1ha[midirdwarfagn]
        o3hb_sel = o3hb[midirdwarfagn]
        ax3.plot(o1ha_sel, o3hb_sel, 'kp', mfc = 'none',
                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
        
        print(100.0*np.sum(dwarfagn)/np.sum(dwarf), \
              100.0*binom_conf_interval(np.sum(dwarfagn),np.sum(dwarf)) - \
              (100.0*np.sum(dwarfagn)/np.sum(dwarf)))
        return (ax1, ax2, ax3)

## -*- coding: utf-8 -*-
#"""
#Created on Tue Mar 22 13:51:04 2022
#
#@author: mugdhapolimera
#"""
#
## -*- coding: utf-8 -*-
#"""
#Created on Fri Mar 26 17:31:14 2021
#
#@author: mugdhapolimera
#"""
#
##This program makes a Line-Ratio diagram (also known as a BPT plot or Kewley diagram)
##with labels using data from the RESOLVE survey to classify galaxies as LINERs,
##Seyferts, Composites, or AGNs on the basis of their flux ratio for distinction.
#
##Original code from Ashley Bittner 08/03/2017
##Edited version from Margie Bruff 01/07/2018
##Updated by Carlynn Ferguson 03/30/2018
#
##suggested use of python debugger to understand the code more thoroughly
##KDE plotting 
##https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density/
##import pdb
#
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import os 
#import sys
#from scipy.stats import kde
#import matplotlib
#import matplotlib.cm as cm
#import matplotlib.colors as colors
#matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#matplotlib.rcParams.update({'font.size': 20})
#matplotlib.rcParams.update({'axes.linewidth': 2})
#matplotlib.rcParams.update({'lines.linewidth': 2})
#from astropy.stats import binom_conf_interval
#
##define demarcation function: log_NII_HA vs. log_OIII_HB
#def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
#    return 1.3 + (0.61 / (log_NII_HA - 0.05))
#def n2halocus(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
#    return 1.1 + (0.61 / (log_NII_HA + 0.08))
#def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
#    return 1.19 + (0.61 / (log_NII_HA - 0.47))
##    return 0.57 + (0.13 / (log_NII_HA - 0.003))
#def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
#    return 1.30 + (0.72 / (log_SII_HA - 0.32))
#    #return 0.58 + (0.04 / (log_SII_HA +0.012))
#def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
#    return 0.76 + 1.89*log_SII_HA
#def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
#    return 1.33 + (0.73 / (log_OI_HA + 0.59))
#    #return 0.61 + (0.056 / (log_OI_HA + 0.40))
#def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
#    return 1.3 + 1.18*log_OI_HA
#def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
#    return -0.59
#
#def he2hbmain(log_NII_HA):
#    return -1.22+1.0/(8.92*log_NII_HA+1.32)
#
#def he2hbmainclass(log_NII_HA):
#    refn2ha = np.linspace(-3.0, -0.15)
#    main = he2hbmain(refn2ha)
#    line = np.poly1d(np.polyfit(refn2ha, main, 15))
#    return line(log_NII_HA) 
#
#def he2hblimit(log_NII_HA):
#        return -1.07+1.0/(8.92*log_NII_HA-0.95)
#
#def he2hblimitclass(log_NII_HA):
#    refn2ha = np.linspace(-3.0, -0.15)
#    limit = he2hblimit(refn2ha)
#    line = np.poly1d(np.polyfit(refn2ha, limit, 15))
#    return line(log_NII_HA) 
#
#def ratioerror(num,num_err,den, den_err):
#    err_num2 = (num_err/(num*np.log(10)))**2
#    err_den2 = (den_err/(den*np.log(10)))**2
#    return np.sqrt(err_num2 + err_den2)
#    # note that the ratio uses only the stronger line, but for S/N reasons we add
#    # the weaker and multiply by 3/4 since Chris Richardson says the canonical
#    # line ratio is 3:1 (this needs to be updated with a more precise number)
#def bpt_plots(inputfile, outputfile,  selinputfile, s06outputfile,  midirfile, 
#              eco, resolve, full, sdsscat, save, ax1, ax2, ax3, simple, allcat):
#    df = pd.read_csv(inputfile)
#        #define alternate catalog names
#    if 'name' in df.keys():
#        df['NAME'] = df['name']
#    if 'CATID' in df.keys():
#        df['NAME'] = df['CATID']
#    name = df['name']
#    df['NAME'] = df['name']
#    df.index = df.name
#
#    s06 = pd.read_csv(s06outputfile)
#    s06.index = s06.galname
##    s06agn = np.array(s06[~s06.defstarform].galname)
#    s06agn = np.array(s06.galname[s06.composite])
#
#    midir = pd.read_csv(midirfile)
#    midir.index = midir.name
#    midiragn = np.array(midir.name.iloc[np.where(midir.agnflag == True)[0]])
#
#    if sdsscat == 'port':
#        nii = df['Flux_NII_6583']
#        nii_sum = (df['Flux_NII_6583']+ df['Flux_NII_6547'])*3./4
#        nii_sum_err = (np.sqrt(df['Flux_NII_6547_Err']**2 + df['Flux_NII_6583_Err']**2))*3./4
#        oiii = df['Flux_OIII_5006']
#        oiii_err = df['Flux_OIII_5006_Err']
#        h_alpha = df['Flux_Ha_6562']
#        h_alpha_err = df['Flux_Ha_6562_Err']
#        h_beta = df['Flux_Hb_4861']
#        h_beta_err = df['Flux_Hb_4861_Err']
#        oi = df['Flux_OI_6300']
#        oi_err = df['Flux_OI_6300_Err']
#        sii_sum = df['Flux_SII_6716'] + df['Flux_SII_6730']
#        sii_sum_err = np.sqrt(df['Flux_SII_6716_Err']**2 + df['Flux_SII_6730_Err']**2)
#    if sdsscat == 'jhu' or sdsscat == 'nsa' or sdsscat == 'master':
#        
#        nii = df['nii_6584_flux']
#        if 'nii_6548_flux' in df.keys():
#            nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
#            nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
#        else:        
#            nii_sum = df['nii_6584_flux']
#            nii_sum_err = df['nii_6584_flux_err']
#
#        oiii = df['oiii_5007_flux']
#        oiii_err = df['oiii_5007_flux_err']
#        h_alpha = df['h_alpha_flux']
#        h_alpha_err = df['h_alpha_flux_err']
#        h_beta = df['h_beta_flux']
#        h_beta_err = df['h_beta_flux_err']
#        oi = df['oi_6300_flux']
#        oi_err = df['oi_6300_flux_err']
#        if 'sii_6717_flux' in df.keys():
#            sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
#    
#            sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
#        else:
#            sii_sum = df['sii_6731_flux']
#        
#            sii_sum_err = df['sii_6731_flux_err']
#    
#    
#    data = df.name > 0
#        
#    if eco: 
#        resname = df['resname'] #for eco
#        resname = resname != 'notinresolve'
#    if resolve:
#        econame = df['NAME']#df['econame'] #for resolve
#        econame = df['NAME']#econame != 'notineco'
#
#    #print total points shared with alternate catalog
#    if full: 
#        sel = (np.where(data)[0]) #for eco
#    elif eco: 
#        sel = (np.where(data & resname)[0]) #for eco
#    elif resolve:
#        sel = (np.where(data & econame)[0]) #for resolve
#    print ''
#    print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)
#    #df = df[data]
#    
#    nii = nii[data]
#    nii_sum = nii_sum[data]
#    nii_sum_err = nii_sum_err[data]
#    oiii = oiii[data]
#    oiii_err = oiii_err[data]
#    oi = oi[data]
#    oi_err = oi_err[data]
#    sii_sum = sii_sum[data]
#    sii_sum_err = sii_sum_err[data]
#    h_beta = h_beta[data]
#    h_beta_err = h_beta_err[data]
#    h_alpha = h_alpha[data]
#    h_alpha_err = h_alpha_err[data]
#    
#    if sdsscat =='nsa':
#        df = df[data]
#    
#        
#    sel = pd.read_csv(selinputfile)
#    sel.index = sel.name
#    onlybpt = set(df.name) - set(sel.name)
#    selname = set(sel.name)
##    oi.loc[onlybpt] = np.nan
##    oi_err.loc[onlybpt] = np.nan
##    sii_sum.loc[onlybpt] = np.nan
##    sii_sum_err.loc[onlybpt] = np.nan
#
#    siisel = pd.read_csv('ECO+RESOLVE_bpt1snr5_siisnr_jhu.csv')
#    siisel.index=siisel.name
#    siisel = siisel['snr']
#    oisel = pd.read_csv('ECO+RESOLVE_bpt1snr5_oisnr_jhu.csv')
#    oisel.index=oisel.name
#    oisel = oisel['snr']
#    
#    #bpt = pd.read_csv('ECO+RESOLVE_bpt1snr5_dext_'+sdsscat+'.csv')
#    #bpt.index = bpt.name
#    #onlybpt = set(bpt.name)
#    #selname = set(df.name) - set(bpt.name)
#    
#
#    
#    #length of data to be used for debugging
#    datalen = np.sum(data)
#    subsetname = df.NAME[data]
#    
#    # data ratios
#    #n2ha = np.log10(nii_sum/h_alpha)
#    o3hb = np.log10(oiii/h_beta) # always the y-axis
#    o1ha = np.log10(oi/h_alpha)
#    s2ha = np.log10(sii_sum/h_alpha)
#    n2ha = np.log10(nii_sum/h_alpha)
#    
#    #Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
#    #and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.
#    
#    #NII plot selectors
#    compsel1 = (o3hb >= n2hacompmin(n2ha)) & (o3hb <= n2hamain(n2ha))
#    sfsel1 = (o3hb < n2hacompmin(n2ha)) & (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
#    agnsel1= (o3hb > n2hamain(n2ha))
#    #plt.hist(o1ha_err[o1ha_err < 1e5], bins = 'fd')
#
#    #SII plot selectors
#    sfsel2 = (o3hb <= s2hamain(s2ha)) & ~compsel1 & siisel
#    seyfsel2 = ((o3hb > s2hamain(s2ha)) & (o3hb >= s2halinseyf(s2ha))) & siisel
#    linersel2 = ((o3hb > s2hamain(s2ha)) & (o3hb < s2halinseyf(s2ha))) & siisel
#    agnsel2 = (o3hb > s2hamain(s2ha)) & ~compsel1 & siisel
#    
#    #OI plot selectors
#    sfsel3 = (o3hb <= o1hamain(o1ha)) & (o1ha < -0.7) & ~compsel1 & oisel
#    seyfsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb >= o1halinseyf(o1ha)) & oisel
#    linersel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb < o1halinseyf(o1ha)) & oisel
#    agnsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & ~compsel1 & oisel
#    
#    #REFERENCE for cumulative plot selectors
#    seyfselr = seyfsel2 & seyfsel3
#    linerselr = linersel2 & linersel3
#    
#    #cumulative plot selectors
#    sfsel = sfsel1 & sfsel2 & sfsel3 
#    compsel = compsel1  #composite galaxies
#    seyfsel = agnsel1 & seyfselr #Seyfert AGN galaxies
#    linersel = agnsel1 & linerselr #LINER AGN galaxies
#    ambigsel1 = sfsel1 & (agnsel2 | agnsel3) & (siisel & oisel) #SF in first plot, AGN in subsequent plot
#    ambigsel2 = np.array(agnsel1) & (np.array(sfsel2) | np.array(sfsel3)) #AGN in first plot, SF in subsequent plot
#    ambagnsel = agnsel1 & ~seyfselr & ~linerselr & ~(sfsel2 | sfsel3) #Ambiguous AGN
#    
#    sftoagn1 = sfsel1 & agnsel2
#    sftoagn2 = sfsel1 & agnsel3
#
#    #Save the BPT flags to a CSV file
#    emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
#    defagn = seyfsel | linersel | ambagnsel
#    convagn = compsel | defagn
#    
#    flags = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
#                              'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
#                              'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
#                              'sftoagn1':sftoagn1, 'sftoagn2': sftoagn2, 'convagn': convagn})
#    if save:
#        flags.to_csv(outputfile ,index=False)
##    keys = ['defagn', 'composite', 'sftoagn', 'agntosf','dwarfagn']
#    keys = ['convagn', 'sftoagn']
#    marker = {'convagn': '1', 'agntosf': 'c^', 'ambigagn': 'ms', 'composite': 'ms', 
#              'defagn': 'ro', 
#              'defliner': 'yo', 'defseyf': 'co', 'dwarfagn': 'ks',
#              'sftoagn': 'bs'}
#    markersize = {'convagn': 10, 'agntosf': 14, 'ambigagn': 8, 'composite': 8, 'defagn': 8, 
#              'defliner': 8, 'defseyf': 8, 'sftoagn': 8, 
#              'dwarfagn':16}
#    alpha = {'convagn': 1, 'agntosf': 1, 'ambigagn': 0.5, 'composite': 0.5, 'defagn': 0.5, 
#              'defliner': 0.5, 'defseyf': 0.5, 'sftoagn': 0.7, 
#              'dwarfagn':1}
#    
#    markercolors = {'convagn': 'r', 'agntosf': 'c', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
#              'defliner': 'y', 'defseyf': 'c', 'dwarfagn': 'k', 'sftoagn': 'b'}
#    
#    labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
#              'composite': 'Composite', 'defagn': 'Traditional AGN', 
#              'defliner': 'LINER', 'defseyf': 'Seyfert', 
#              'dwarfagn': 'SEL Dwarf AGN', 'defstarform': 'Definite SF', 
#              'sftoagn': 'SF-AGN', 'convagn': 'Conventional AGN'}
#    
#    
#    #checking that plotted points are within the total data range
#    print ''
#    sfselpts = (len(np.where(sfsel)[0]))
#    seyfselpts = (len(np.where(seyfsel)[0]))
#    linerselpts = (len(np.where(linersel)[0]))
#    compselpts = (len(np.where(compsel)[0]))
#    agnselpts = (len(np.where(ambagnsel)[0]))
#    ambigsel1pts = (len(np.where(ambigsel1)[0]))
#    ambigsel2pts = (len(np.where(ambigsel2)[0]))
#    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+\
#    ambigsel1pts+ambigsel2pts
#    sfpercent = float(sfselpts)/float(datalen)*100
#    seyfpercent = float(seyfselpts)/float(datalen)*100
#    linerpercent = float(linerselpts)/float(datalen)*100
#    comppercent = float(compselpts)/float(datalen)*100
#    agnpercent = float(agnselpts)/float(datalen)*100
#    ambig1percent = float(ambigsel1pts)/float(datalen)*100
#    ambig2percent = float(ambigsel2pts)/float(datalen)*100
#    #df.index = fulldf.name
#    if 'logmstar' in df.keys():
#        dwarf = (df.logmstar < 9.5)
#        giant = (df.logmstar > 9.5)
#    
#    agn = (ambigsel1|seyfsel|linersel|ambagnsel|compsel|ambigsel2)
#    dwarfagn = dwarf & agn
#    giantagn = giant & agn
#
#    s06dwarfagn = np.intersect1d(s06agn, df.name[dwarf])
#    midirdwarfagn = np.intersect1d(midiragn, df.name[dwarf])
#    
#    print ("DATA POINTS: "),datalen
#    print ("TOTAL PLOTTED POINTS: "), totalselpts
#    print ("TOTAL PLOTTED POINTS OMITTED: "), datalen-totalselpts
#    print "* IF NUMBERS ABOVE ARE AT ALL NEGATIVE THEN THERE IS OVERPLOTTING"
#    print ("Definite Star Forming: "),sfselpts,("("),round(sfpercent, 2),("%"),(")")
#    print ("Composite: "),compselpts, ("("),round(comppercent, 2),("%"),(")")
#    print ("SF --> AGN: "), ambigsel1pts, ("("),round(ambig1percent, 2),("%"),(")")
#    print ("AGN --> SF: "), ambigsel2pts, ("("),round(ambig2percent, 2),("%"),(")")
#    print ("Ambiguous AGN: "),agnselpts, ("("),round(agnpercent, 2),("%"),(")")
#    print ("Seyfert: "),seyfselpts, ("("),round(seyfpercent, 2),("%"),(")")
#    print ("LINER: "),linerselpts, ("("),round(linerpercent, 2),("%"),(")")
#    print ("TOTAL KNOWN AGN: "),linerselpts+seyfselpts+agnselpts, ("("), \
#    round(linerpercent+seyfpercent+agnpercent, 2), ("% )")
#    print ("POSSIBLE TOTAL AGN: "),linerselpts+seyfselpts+agnselpts+ambigsel1pts+ambigsel2pts,("("),\
#    round(linerpercent+seyfpercent+agnpercent+ambig1percent+ambig2percent, 2), ("% )")
#    print ("Percent Omitted: "), round((100-(sfpercent+seyfpercent+linerpercent+comppercent+agnpercent+ambig1percent+ambig2percent)), 2), ("%")
#    print ''
#    
#    print ("AGN in Dwarf Galaxies: "), 100*round(np.sum(dwarfagn)/float(np.sum(dwarf)),2), ("%")
#    print ("AGN in Giant Galaxies: "), 100*round(np.sum(giantagn)/float(np.sum(giant)),2), ("%")
#    print ("AGN in dwarfs: "), np.sum(agn & dwarf)
#    print ("Number of Dwarfs:"), np.sum(dwarf)
#    
#    ###PLOTS###
#    #reference points in x-direction for demarcation lines on plots
#    if simple:
#        refn2ha = np.linspace(-3.0, 0.35)
#        refoiha = np.linspace(-2.5, -0.4)
#        refsiiha = np.linspace(-2, 0.3,100)
#        
#        def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
#          	new_cmap = colors.LinearSegmentedColormap.from_list(
#                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#                cmap(np.linspace(minval, maxval, n)))
#          	return new_cmap
#        sf_colors_map = truncate_colormap(cm.gray_r)
#        ndx = []#'ECO03494']#np.where((df.NAME == 'rs1105') | (df.NAME == 'rs1375'))[0]
#        xmin = refn2ha.min(); xmax = refn2ha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        
#        #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
#        
#        #NII/OIII plot
#        definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#        ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
#
#        ax1.set_xlim(-1.5,0.5)
#        ax1.set_ylim(-1.0,1.0)
#        ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
#        ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#        
#        for key in keys:
#            n2ha_sel = n2ha.loc[flags[key]]
#            o3hb_sel = o3hb.loc[flags[key]]
#            ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
#                     alpha = alpha[key], markersize = markersize[key], 
#                     mew = 2, mec = markercolors[key], label = labels[key])
#
#        if allcat:
#            n2ha_sel = n2ha[dwarfagn]
#            o3hb_sel = o3hb[dwarfagn]
#            ax1.plot(n2ha_sel, o3hb_sel, 'ks', mfc = 'none',
#                     markersize = 14, mew = 2, label = 'Dwarf AGN')
#    
#            n2ha_sel = n2ha[s06dwarfagn]
#            o3hb_sel = o3hb[s06dwarfagn]
#            ax1.plot(n2ha_sel, o3hb_sel, 'v',color = 'lime',  mfc = 'none',
#                     markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')
#    
#            n2ha_sel = n2ha[midirdwarfagn]
#            o3hb_sel = o3hb[midirdwarfagn]
#            ax1.plot(n2ha_sel, o3hb_sel, 'p', color = 'orange', zorder = 10, mfc = 'none',
#                     markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        ax1.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
#        
#        main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Theoretical Maximum Starburst Line (Ke01)')
#        composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
#        'k--', label = 'Composite Line (Ka03)')
#        
#        #SII plot
#        xmin = refsiiha.min(); xmax = refsiiha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        s2ha = s2ha.loc[siisel]
#        o3hb_orig = o3hb.copy()
#        o3hb = o3hb.loc[siisel]
#        
#        s2ha = s2ha[~np.isnan(s2ha)]
#        o3hb = o3hb[~np.isnan(s2ha)]
#        
#        definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#        agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
#        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
#        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
#        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
#                                        ymin_agn:ymax_agn:nbins*1j]
#        k = kde.gaussian_kde(agn_contour.T)
#        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#        ax2.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map)
#        main1, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
#        ax2.set_xlim(-1.5,0.5)
#        ax2.set_ylim(-1,1)
#        ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#        ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
#                          'k-.', label = 'Liner/Seyfert Division')
#        for key in keys:
#            s2ha_sel = s2ha.loc[flags[key]]
#            o3hb_sel = o3hb.loc[flags[key]]
#            ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
#                     alpha = alpha[key], markersize = markersize[key], 
#                     mew = 2, mec = markercolors[key], label = labels[key])
#        s2ha_sel = s2ha[sftoagn1 & ~oisel & seyfsel2]
#        o3hb_sel = o3hb[sftoagn1 & ~oisel & seyfsel2]
#        ax2.plot(s2ha_sel, o3hb_sel, 'ko', mfc = 'none',
#                 markersize = 14, mew = 2, label = 'Dwarf AGN')
#        if allcat:
#            s2ha_sel = s2ha[dwarfagn]
#            o3hb_sel = o3hb[dwarfagn]
#            ax2.plot(s2ha_sel, o3hb_sel, 'ks', mfc = 'none',
#                     markersize = 14, mew = 2, label = 'Dwarf AGN')
#    
#            s2ha_sel = s2ha[s06dwarfagn]
#            o3hb_sel = o3hb[s06dwarfagn]
#            ax2.plot(s2ha_sel, o3hb_sel, 'v', color = 'lime', mfc = 'none',
#                     markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')
#    
#            s2ha_sel = s2ha[midirdwarfagn]
#            o3hb_sel = o3hb[midirdwarfagn]
#            ax2.plot(s2ha_sel, o3hb_sel, 'p', mfc = 'none', color = 'orange', zorder = 10, 
#                     markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        #OI Plot
#        xmin = refoiha.min(); xmax = 0#refoiha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        
#        o1ha = o1ha.loc[oisel]
#        o3hb = o3hb_orig.loc[oisel]
#
#        o1ha = o1ha[~np.isnan(o1ha) & np.isfinite(o1ha)]
#        o3hb = o3hb[~np.isnan(o1ha) & np.isfinite(o1ha)]
#        
#        definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#
#        ax3.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map)
#        main1, = ax3.plot(refoiha[refoiha < -0.75], o1hamain(refoiha[refoiha < -0.75]), 'k', 
#                          label = 'Ke01 Maximum Starburst Line')
#        ax3.set_xlim(-2.0,-0.4)
#        ax3.set_ylim(-1,1)
#        ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
#                                       'k-.', label = 'Ke06 Liner/Seyfert Division Line')
#        ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#
#        for key in keys:
#            o1ha_sel = o1ha.loc[flags[key]]
#            o3hb_sel = o3hb.loc[flags[key]]
#            ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], 
#                     alpha = alpha[key], markersize = markersize[key], 
#                     mew = 2, mec = markercolors[key], label = labels[key])
#        o1ha_sel = o1ha[sftoagn2 & ~siisel]
#        o3hb_sel = o3hb[sftoagn2 & ~siisel]
#        ax3.plot(o1ha_sel, o3hb_sel, 'gs', mfc = 'g',
#                 markersize = 14, mew = 2, label = 'Dwarf AGN')
#        if allcat:
#            o1ha_sel = o1ha[dwarfagn]
#            o3hb_sel = o3hb[dwarfagn]
#            ax3.plot(o1ha_sel, o3hb_sel, 'ks', mfc = 'none',
#                     markersize = 14, mew = 2, label = 'Dwarf AGN')
#    
#            o1ha_sel = o1ha[s06dwarfagn]
#            o3hb_sel = o3hb[s06dwarfagn]
#            ax3.plot(o1ha_sel, o3hb_sel, 'v', color = 'lime', mfc = 'none',
#                     markersize = 14, mew = 2, label = 'S06 Bonus Dwarf AGN')
#    
#            o1ha_sel = o1ha[midirdwarfagn]
#            o3hb_sel = o3hb[midirdwarfagn]
#            ax3.plot(o1ha_sel, o3hb_sel, 'p', mfc = 'none', color = 'orange', zorder = 10, 
#                     markersize = 14, mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        print(100.0*np.sum(dwarfagn)/np.sum(dwarf), \
#              100.0*binom_conf_interval(np.sum(dwarfagn),np.sum(dwarf)) - \
#              (100.0*np.sum(dwarfagn)/np.sum(dwarf)))
#        return (ax1, ax2, ax3)
#        
#    else:
#        refn2ha = np.linspace(-3.0, 0.35)
#        refoiha = np.linspace(-2.5, -0.4)
#        refsiiha = np.linspace(-2, 0.3,100)
#        
#        def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
#          	new_cmap = colors.LinearSegmentedColormap.from_list(
#                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#                cmap(np.linspace(minval, maxval, n)))
#          	return new_cmap
#        sf_colors_map = truncate_colormap(cm.gray_r)
#        ndx = []#'ECO03494']#np.where((df.NAME == 'rs1105') | (df.NAME == 'rs1375'))[0]
#        xmin = refn2ha.min(); xmax = refn2ha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        
#        #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
#        
#        #NII/OIII plot
#        definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#
#        ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
#        ax1.set_xlim(-1.5,0.5)
#        ax1.set_ylim(-1.0,1.0)
#        ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
#        ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#        
#        for key in keys:
#            if key != 'dwarfagn':
#                n2ha_sel = n2ha.loc[flags[key]]
#                o3hb_sel = o3hb.loc[flags[key]]
#                ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
#                         markersize = markersize[key], mew = 0, label = labels[key])
#            else:
#                n2ha_sel = n2ha[dwarfagn]
#                o3hb_sel = o3hb[dwarfagn]
#                ax1.plot(n2ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
#                         markersize = markersize[key], mew = 2, label = labels[key])
#        n2ha_sel = n2ha[s06dwarfagn]
#        o3hb_sel = o3hb[s06dwarfagn]
#        ax1.plot(n2ha_sel, o3hb_sel, 'kv', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
#        n2ha_sel = n2ha[midirdwarfagn]
#        o3hb_sel = o3hb[midirdwarfagn]
#        ax1.plot(n2ha_sel, o3hb_sel, 'kp', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        ax1.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
#        
#        main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Theoretical Maximum Starburst Line (Ke01)')
#        composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
#        'k--', label = 'Composite Line (Ka03)')
#        
#        #SII plot
#        xmin = refsiiha.min(); xmax = refsiiha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        
#        definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#        agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
#        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
#        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
#        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
#                                        ymin_agn:ymax_agn:nbins*1j]
#        k = kde.gaussian_kde(agn_contour.T)
#        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#        ax2.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map)
#        main1, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
#        ax2.set_xlim(-1.5,0.5)
#        ax2.set_ylim(-1,1)
#        ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#        ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
#                          'k-.', label = 'Liner/Seyfert Division')
#        for key in keys:
#            if key != 'dwarfagn':
#                s2ha_sel = s2ha.loc[flags[key]]
#                o3hb_sel = o3hb.loc[flags[key]]
#                ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
#                         markersize = markersize[key], mew = 0, label = labels[key])
#            else:
#                s2ha_sel = s2ha[dwarfagn]
#                o3hb_sel = o3hb[dwarfagn]
#                ax2.plot(s2ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
#                         markersize = markersize[key], mew = 2, label = labels[key])
#        s2ha_sel = s2ha[s06dwarfagn]
#        o3hb_sel = o3hb[s06dwarfagn]
#        ax2.plot(s2ha_sel, o3hb_sel, 'kv', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
#        s2ha_sel = s2ha[midirdwarfagn]
#        o3hb_sel = o3hb[midirdwarfagn]
#        ax2.plot(s2ha_sel, o3hb_sel, 'kp', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        #OI Plot
#        xmin = refoiha.min(); xmax = 0#refoiha.max()
#        ymin = -1.25; ymax = 1.5
#        nbins = 50
#        
#        definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
#        xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#        k2 = kde.gaussian_kde(definite.T)
#        definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#        agn_contour = np.column_stack((o1ha[defagn], o3hb[defagn]))
#        xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
#        ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
#        xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
#                                        ymin_agn:ymax_agn:nbins*1j]
#        k = kde.gaussian_kde(agn_contour.T)
#        agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
#        ax3.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=sf_colors_map)
#        main1, = ax3.plot(refoiha[refoiha < -0.75], o1hamain(refoiha[refoiha < -0.75]), 'k', 
#                          label = 'Ke01 Maximum Starburst Line')
#        ax3.set_xlim(-2.0,-0.4)
#        ax3.set_ylim(-1,1)
#        ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
#                                       'k-.', label = 'Ke06 Liner/Seyfert Division Line')
#        ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#        for key in keys:
#            if key != 'dwarfagn':
#                o1ha_sel = o1ha.loc[flags[key]]
#                o3hb_sel = o3hb.loc[flags[key]]
#                ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], alpha = alpha[key], 
#                         markersize = markersize[key], mew = 0, label = labels[key])
#            else:
#                o1ha_sel = o1ha[dwarfagn]
#                o3hb_sel = o3hb[dwarfagn]
#                ax3.plot(o1ha_sel, o3hb_sel, marker[key], color = markercolors[key], mfc = 'none',
#                         markersize = markersize[key], mew = 2, label = labels[key])
#        o1ha_sel = o1ha[s06dwarfagn]
#        o3hb_sel = o3hb[s06dwarfagn]
#        ax3.plot(o1ha_sel, o3hb_sel, 'kv', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'S06-BPT Dwarf AGN')
#        o1ha_sel = o1ha[midirdwarfagn]
#        o3hb_sel = o3hb[midirdwarfagn]
#        ax3.plot(o1ha_sel, o3hb_sel, 'kp', mfc = 'none',
#                 markersize = markersize[key], mew = 2, label = 'Mid-IR Dwarf AGN')
#        
#        print(100.0*np.sum(dwarfagn)/np.sum(dwarf), \
#              100.0*binom_conf_interval(np.sum(dwarfagn),np.sum(dwarf)) - \
#              (100.0*np.sum(dwarfagn)/np.sum(dwarf)))
#        return (ax1, ax2, ax3)
#>>>>>>> 4d300bdd936fd05600c51bbe79ca4a5a518e4724
