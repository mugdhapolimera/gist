# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 11:37:13 2021

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May 26 15:38:42 2019

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
matplotlib.rcParams.update({'font.size': 15})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import binom_conf_interval

#Read in RESOLVE/ECO extinction corrected and S/N filtered data
eco = 0
resolve = 0
save = 0
full = 1
#sdsscat = 'jhu'

#sdsscat = 'port'
sdsscat = 'nsa'
#sdsscat = 'master'
if sys.platform == 'linux2':
        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')

else:
    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')


inputfile = 'ECO+RESOLVE_snr5_dext_'+sdsscat+'.csv'
print 'ECO+RESOLVE RESULTS'
outputfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
df = pd.read_csv(inputfile)
if 'source' not in df.keys():
    df['source'] = sdsscat

#define alternate catalog names
if 'name' in df.keys():
    df['NAME'] = df['name']
#if (sdsscat == 'nsa') & ('resname' in df.keys()):
#    df['name'] = df['resname']
if 'CATID' in df.keys():
    df['NAME'] = df['CATID']
name = df['name']
df['NAME'] = df['name']
df.index = df.name
#df = df.loc[econsanames]
#df = df.loc[inobssample.index.values[inobssample]]

if he2_flag:
    df = df[~np.isnan(df.Flux_HeII_4685)]
    print(len(df))
    heii = df['Flux_HeII_4685']
    heii_err = df['Flux_HeII_4685_Err']
#df['name'] = df['resname']
if eco: 
    resname = df['resname'] #for eco
    resname = resname != 'notinresolve'
if resolve:
    econame = df['NAME']#df['econame'] #for resolve
    econame = df['NAME']#econame != 'notineco'

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
    heii = df['Flux_HeII_4685']
    heii_err = df['Flux_HeII_4685_Err']
if sdsscat == 'jhu' or sdsscat == 'nsa' or sdsscat == 'master':
    
    nii = df['nii_6584_flux']
    if 'nii_6548_flux' in df.keys():
        nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
        nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
    else:        
        nii_sum = df['nii_6584_flux']
        nii_sum_err = df['nii_6584_flux_err']
    #nii_sum[df.source == 'nsa'] = df['nii_6584_flux']
    #nii_sum_err[df.source == 'nsa'] = df['nii_6584_flux_err']

    #nii = nii_sum
    # note that the ratio uses only the stronger line, but for S/N reasons we add
    # the weaker and multiply by 3/4 since Chris Richardson says the canonical
    # line ratio is 3:1 (this needs to be updated with a more precise number)
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
#    sii_sum[df.source == 'nsa'] = df['sii_6731_flux']
#    sii_sum_err[df.source == 'nsa'] = df['sii_6731_flux_err']

#nii = df['Flux_NII_6583']
#nii_sum = (df['Flux_NII_6583']+ df['Flux_NII_6547'])*3./4
#nii_sum_err = (np.sqrt(df['Flux_NII_6547_Err']**2 + df['Flux_NII_6583_Err']**2))*3./4
#oiii = df['Flux_OIII_5006']
#oiii_err = df['Flux_OIII_5006_Err']
#h_alpha = df['Flux_Ha_6562']
#h_alpha_err = df['Flux_Ha_6562_Err']
#h_beta = df['Flux_Hb_4861']
#h_beta_err = df['Flux_Hb_4861_Err']
#oi = df['Flux_OI_6300']
#oi_err = df['Flux_OI_6300_Err']
#sii_sum = df['Flux_SII_6716'] + df['Flux_SII_6730']
#sii_sum_err = np.sqrt(df['Flux_SII_6716_Err']**2 + df['Flux_SII_6730_Err']**2)

#Filter Data: all non-negative SEL fluxes and errors; Hbeta >3sigma
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
            (sii_sum > 0) & (h_beta > 0) & (h_beta > 5*h_beta_err) &
            (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0) & 
            (oi_err > 0) & (sii_sum_err > 0))

snr = ((h_alpha > 5*h_alpha_err) & (nii_sum > 5*nii_sum_err) & (oiii > 5*oiii_err) & 
       (oi > 5*oi_err) & (sii_sum > 5*sii_sum_err) & (h_beta > 5*h_beta_err))

#if he2_flag:
#    he2data = (heii/heii_err >=3) & (heii_err > 0)
#    gooddata = (h_beta > 0) & (h_beta > 5*h_beta_err) & \
#                (h_alpha_err > 0) & (nii_sum_err > 0) & \
#                (h_alpha > 0) & (nii_sum > 0)
#    snr = (nii_sum> 5*nii_sum_err) & (h_alpha > 5*h_alpha_err)
#    data = he2data
#    print(np.sum(he2data),np.sum(data))
#else:
data = df.name > 0 #gooddata #& snr #use ALL galaxy data within catalog
if sdsscat == 'master':
    data = df.name > 0
#print total points shared with alternate catalog
if full: 
    sel = (np.where(data)[0]) #for eco
elif eco: 
    sel = (np.where(data & resname)[0]) #for eco
elif resolve:
    sel = (np.where(data & econame)[0]) #for resolve
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

#nii = df['nii_6584_flux']
##if 'nii_6548_flux' in df.keys():
#nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
#nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
##else:
## note that the ratio uses only the stronger line, but for S/N reasons we add
## the weaker and multiply by 3/4 since Chris Richardson says the canonical
## line ratio is 3:1 (this needs to be updated with a more precise number)
#oiii = df['oiii_5007_flux']
#oiii_err = df['oiii_5007_flux_err']
#h_alpha = df['h_alpha_flux']
#h_alpha_err = df['h_alpha_flux_err']
#h_beta = df['h_beta_flux']
#h_beta_err = df['h_beta_flux_err']
#oi = df['oi_6300_flux']
#oi_err = df['oi_6300_flux_err']
#if 'sii_6717_flux' in df.keys():
#    sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
#
#    sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
#else:
#    sii_sum = df['sii_6731_flux']
#
#    sii_sum_err = df['sii_6731_flux_err']

if he2_flag:
    heii = heii[data] # 3-sigma cut for HeII selection
subsetname = df.NAME[data]
if sdsscat =='nsa':
    df = df[data]
#    df.to_csv('NSA_ECO_snr5.csv')
#length of data to be used for debugging
datalen = np.sum(data)

# data ratios
#n2ha = np.log10(nii_sum/h_alpha)
o3hb = np.log10(oiii/h_beta) # always the y-axis
o1ha = np.log10(oi/h_alpha)
s2ha = np.log10(sii_sum/h_alpha)
n2ha = np.log10(nii_sum/h_alpha)
s2ha_err = np.array(ratioerror(sii_sum, sii_sum_err, h_alpha, h_alpha_err))
o3hb_err = np.array(ratioerror(oiii, oiii_err, h_beta, h_beta_err))
o1ha_err = np.array(ratioerror(oi, oi_err, h_alpha, h_alpha_err))
n2ha_err = np.array(ratioerror(nii_sum, nii_sum_err, h_alpha, h_alpha_err))

if he2_flag:
    he2hb = np.log10(heii/h_beta)

#Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
#and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.

#NII plot selectors
compsel1 = (o3hb >= n2hacompmin(n2ha)) & (o3hb <= n2hamain(n2ha))
sfsel1 = (o3hb < n2hacompmin(n2ha)) & (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
agnsel1= (o3hb > n2hamain(n2ha))
#plt.hist(o1ha_err[o1ha_err < 1e5], bins = 'fd')
#SII plot selectors
sfsel2 = (o3hb <= s2hamain(s2ha)) & ~compsel1
seyfsel2 = ((o3hb > s2hamain(s2ha)) & (o3hb >= s2halinseyf(s2ha)))
linersel2 = ((o3hb > s2hamain(s2ha)) & (o3hb < s2halinseyf(s2ha)))
agnsel2 = (o3hb > s2hamain(s2ha)) & ~compsel1

#OI plot selectors
sfsel3 = (o3hb <= o1hamain(o1ha)) & (o1ha < -0.7) & ~compsel1
seyfsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb >= o1halinseyf(o1ha))
linersel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb < o1halinseyf(o1ha))
agnsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & ~compsel1

#HEII plot selectors
if he2_flag:
    refn2ha = np.linspace(-3.0, 0.35)
    main = he2hbmain(refn2ha[refn2ha < -0.15])
    limit = he2hblimit(refn2ha[refn2ha < -0.15])
#    sfsel4 = (he2hb <= main) & (he2hb <= limit)
#    agnsel4 = (he2hb > main) & (he2hb > limit)
    sfsel4 = (he2hb <= he2hbmain(n2ha)) & \
                (he2hb <= he2hblimit(n2ha))
    agnsel4 = ((he2hb > he2hbmain(n2ha)) | (n2ha >= -0.2)) & \
                ((he2hb > he2hblimit(n2ha)) | (n2ha >= 0))
#for BPT comparison

#REFERENCE for cumulative plot selectors
seyfselr = seyfsel2 & seyfsel3
linerselr = linersel2 & linersel3

#cumulative plot selectors
if he2_flag:
    sfsel = sfsel1 & sfsel2 & sfsel3 & sfsel4 #definite star forming galaxies
else:
    sfsel = sfsel1 & sfsel2 & sfsel3 
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
if not he2_flag:    
    flags = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
                          'sftoagn1':sftoagn1, 'sftoagn2': sftoagn2})
else:
    flags = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
                          'heiisel':agnsel4})
        
if save:
    flags.to_csv(outputfile ,index=False)

#checking that plotted points are within the total data range
print ''
sfselpts = (len(np.where(sfsel)[0]))
seyfselpts = (len(np.where(seyfsel)[0]))
linerselpts = (len(np.where(linersel)[0]))
compselpts = (len(np.where(compsel)[0]))
agnselpts = (len(np.where(ambagnsel)[0]))
ambigsel1pts = (len(np.where(ambigsel1)[0]))
ambigsel2pts = (len(np.where(ambigsel2)[0]))
if he2_flag:
    heiiselpts = (len(np.where(agnsel4)[0]))
    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+\
    ambigsel1pts+ambigsel2pts+heiiselpts
else:
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

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharey=True)
ax4.set_axis_off()

definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((n2ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
ax1.set_xlim(-1.5,0.5)
ax1.set_ylim(-1.0,1.0)
#sfdata1, = ax1.plot(n2ha[sfsel], o3hb[sfsel], '.', color = 'gray', 
#                        markersize = 8, mew = 0, label = 'SF')
ambig1data1, = ax1.plot(n2ha[ambigsel1], o3hb[ambigsel1], 'bs',
                alpha = 0.5, markersize = 8, mew = 0, label = 'SF-AGN')
compdata1, = ax1.plot(n2ha[compsel], o3hb[compsel], 'ms', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
defagndata1, = ax1.plot(n2ha[defagn], o3hb[defagn], 'rs', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Definite AGN')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3, 
#            colors='k', corner_mask = True, extend = 'both')

if he2_flag:
    agndata4, = ax1.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8,
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
ambig2data1, = ax1.plot(n2ha[ambigsel2], o3hb[ambigsel2],'c^', markersize = 10,
                        mew = 1, label = 'Low-[S II] AGN')
ambig1data1, = ax1.plot(n2ha[dwarf&ambigsel1], o3hb[dwarf&ambigsel1], 'bs',
                markersize = 14, mew = 2, mec = 'k', 
                mfc = 'none', label = 'Dwarf AGN')
#ax1.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Theoretical Maximum Starburst Line (Ke01)')
composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
'k--', label = 'Composite Line (Ka03)')
#ambig1data1, = ax1.plot(n2ha[dwarf&ambigsel1], o3hb[dwarf&ambigsel1], 'bs',
#                markersize = 14, mew = 2, mec = 'k', 
#                mfc = 'none', label = 'Dwarf AGN')
compdata1, = ax1.plot(n2ha[dwarf&compsel], o3hb[dwarf&compsel], 'ms', 
                      markersize = 14, mew = 2, mec = 'k', 
                      mfc = 'none', label = 'Composite')
defagndata1, = ax1.plot(n2ha[dwarf&defagn], o3hb[dwarf&defagn], 'rs', 
                      markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'Definite AGN')
ambig2data1, = ax1.plot(n2ha[dwarf&ambigsel2], o3hb[dwarf&ambigsel2],'cs', 
                        markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'AGN-to-SF')
ambig2data1, = ax1.plot(n2ha[ambigsel2], o3hb[ambigsel2],'c^', markersize = 10,
                        mew = 1, label = 'Low-[S II] AGN')

xmin = refsiiha.min(); xmax = refsiiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50

#fig = plt.figure()
#ax2 = fig.add_subplot(132, sharey = True)

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
#ax2.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3,
#            colors='k')
ax2.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=sf_colors_map)
ax2.plot(s2ha[ambigsel1], o3hb[ambigsel1],
                            'bs', alpha = 0.5, 
                            ms = 8, mew = 0, label = 'SF-AGN')
compdata1, = ax2.plot(s2ha[compsel], o3hb[compsel], 'ms', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
agn1, = ax2.plot(s2ha[defagn], o3hb[defagn], 'rs', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Traditional AGN')
ambig2data1, = ax2.plot(s2ha[ambigsel2], o3hb[ambigsel2],'c^', markersize = 10,
                        mew = 1, label = 'Low-[S II] AGN')
ambig1data1, = ax2.plot(s2ha[dwarf&ambigsel1], o3hb[dwarf&ambigsel1], 'bs',
                markersize = 14, mew = 2, mec = 'k', 
                mfc = 'none', label = 'Dwarf AGN (any type)')
ax2.legend(loc=3, numpoints = 1, fontsize = 15)#, fontsize = 14)
main1, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
ax2.set_xlim(-1.5,0.5)
ax2.set_ylim(-1,1)
compdata1, = ax2.plot(s2ha[dwarf&compsel], o3hb[dwarf&compsel], 'ms', 
                      markersize = 14, mew = 2, mec = 'k', 
                      mfc = 'none', label = 'Composite')
defagndata1, = ax2.plot(s2ha[dwarf&defagn], o3hb[dwarf&defagn], 'rs', 
                      markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'Definite AGN')
ambig2data1, = ax2.plot(s2ha[dwarf&ambigsel2], o3hb[dwarf&ambigsel2],'c^', 
                        markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'AGN-to-SF')
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                  'k-.', label = 'Liner/Seyfert Division')

if he2_flag:
    agndata4, = ax2.plot(s2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, 
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')

#OI Plot
xmin = refoiha.min(); xmax = 0#refoiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50
#fig = plt.figure()
#ax3 = fig.add_subplot(133, sharey = True)

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
#ax3.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3, 
#            colors='k')
ax3.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=sf_colors_map)
composite_fake, = ax3.plot(refoiha[refoiha < -0.75], 
                           o1hamain(refoiha)[refoiha < -0.75], 'k--', 
                           label = 'Ka03 Composite Line')
main1, = ax3.plot(refoiha[refoiha < -0.75], o1hamain(refoiha[refoiha < -0.75]), 'k', 
                  label = 'Ke01 Maximum Starburst Line')
ax3.set_xlim(-2.0,-0.4)
ax3.set_ylim(-1,1)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
compdata1, = ax3.plot(o1ha[compsel], o3hb[compsel], 'ms', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
defagndata1, = ax3.plot(o1ha[defagn], o3hb[defagn], 'rs', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Definite AGN')
#ax3.errorbar(o1ha[ambigsel1], o3hb[ambigsel1], xerr = o1ha_err[ambigsel1],
#                            yerr = o3hb_err[ambigsel1], fmt = 'b.', alpha = 0.5,
#                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
ax3.plot(o1ha[ambigsel1], o3hb[ambigsel1], 'bs', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                               'k-.', label = 'Ke06 Liner/Seyfert Division Line')
if he2_flag:
    agndata4, = ax3.plot(o1ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, 
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
ambig1data1, = ax3.plot(o1ha[dwarf&ambigsel1], o3hb[dwarf&ambigsel1], 'bs',
                markersize = 14, mew = 2, mec = 'k', 
                mfc = 'none', label = 'SF-AGN')
compdata1, = ax3.plot(o1ha[dwarf&compsel], o3hb[dwarf&compsel], 'ms', 
                      markersize = 14, mew = 2, mec = 'k', 
                      mfc = 'none', label = 'Composite')
defagndata1, = ax3.plot(o1ha[dwarf&defagn], o3hb[dwarf&defagn], 'rs', 
                      markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'Definite AGN')
ambig2data1, = ax3.plot(o1ha[dwarf&ambigsel2], o3hb[dwarf&ambigsel2],'c^', 
                        markersize = 14, mew = 2, mec = 'k', mfc = 'none', label = 'Low-[S II] AGN')
ambig2data1, = ax3.plot(o1ha[ambigsel2], o3hb[ambigsel2],'c^', markersize = 12,
                        mew = 1,  label = 'AGN-to-SF')

sftoagn = df[ambigsel1 & dwarf][['radeg','dedeg']]
c = SkyCoord(ra = list(sftoagn.radeg)*u.degree, 
             dec=list(sftoagn.dedeg)*u.degree, frame='icrs')
sftoagn['h'],sftoagn['m'],sftoagn['s'] = c.ra.hms
#print(sftoagn)
#sftoagn.to_csv(sdsscat+'_sfagn.csv')
print(100.0*np.sum(dwarfagn)/np.sum(dwarf), \
      100.0*binom_conf_interval(np.sum(dwarfagn),np.sum(dwarf)) - \
      (100.0*np.sum(dwarfagn)/np.sum(dwarf)))

plt.show()
