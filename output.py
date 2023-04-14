# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 18:05:54 2019

@author: mugdhapolimera
"""

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.table import Table
import os
from pylab import *
import pandas as pd
from scipy.io import readsav
import math

#folder = 'F:/mugdhapolimera/Shared VM/binned3drf0477crop_snr3test2/'binned3drf0477crop_snr10/'
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59
def ratioerror(num,num_err,den, den_err):
    err = (num/den) * np.sqrt((num_err/num)**2 + (den_err/den)**2)
    return err

#galname = 'rf0376'
galname = 'rs0183' #strong OI
#galname = 'rf0503' #strong OI - very compact; only bins 0 and 2 have good SNR
#galname = 'rf0477' #- strong OI; good SNR ~ 7
#galname = 'rs0105' - strong OI but weak SNR for spectrum
#galname = 'rs1143'
#galname = "rf0078"
#galname = 'rs0124' #very low SNR throughout : GIST coulg not bin upto SNR ~ 5

#inputname = "rf0078_finalcube"
inputname = "binned3d"+galname+"crop"

#Directory where GIST is run from
root1 = '/srv/two/sheila/mugpol/github/gist_wrappers/rungist/'
data = fits.open(root1+"inputData/"+inputname+'.fits')[0].data
hgal = fits.open(root1+"inputData/"+inputname+'.fits')[0].header

#Directory where results are stored
root = '/srv/two/sheila/mugpol/github/gist_wrappers/rungist/results/'
folder = galname+'_snr30newdeg'  #name of the results folder
os.chdir(root+folder)

#Continuum fit from ppxf
cfit = Table.read(folder+'_kin-bestfit.fits')
#Try to read continuum+emission line fit fro gandalf if it exists
try: 
    fit = Table.read(folder+'_gas-bestfit_BIN.fits')
except IOError:
    print("No Emission line module")

#Wavelength solution of the input spectrum
lamdata = (np.arange(hgal['NAXIS3']) + hgal['CRPIX3']-1) * hgal['CDELT3'] + hgal['CRVAL3'] 
#Table with info about bins and SNR
table = fits.open(folder+'_table.fits')[1].data
#Wavelength of the binned spectrum from GIST
lam = fits.open(folder+'_BinSpectra.fits')[2].data.LOGLAM
#Flux/intensity of the binned spectrum from GIST
spec = fits.open(folder+'_BinSpectra.fits')[1].data.SPEC

#Wavelength array of the best fit continuum template from GIST
plam = fits.open(folder+'_kin-optimalTemplates.fits')[2].data['LOGLAM_TEMPLATE']
#Wavelength region that is fit (can be read in from config file)
good = (lamdata > 4300) & (lamdata < 6900)
#Best fit continuum template
temp = fits.open(folder+'_kin-optimalTemplates.fits')

#Table with kinematic data
kin = fits.open(folder+'_kin.fits')[1].data
v_gist = kin['V']
sig_gist = kin["SIGMA"][0]
sig_gist_err = kin["FORM_ERR_SIGMA"][0]

#Bin numbers given in the table
nbins = np.unique(np.abs(table.BIN_ID)) #[0]
#Read in RESOLVE public database
df= pd.read_csv('~/github/SDSS_spectra/RESOLVE_filter_new.csv')
df.index = df.name    
#Read in RESOLVE internal database
db = readsav('/srv/one/resolve/database_internal/merged_idl_catalog/stable/resolvecatalog.dat')

#Calculate the radius of the galaxy in pc 
name = [x.decode('utf-8') for x in db['name']]
vhel = dict(zip(name,db['vhel'])) #km/s
v = vhel[galname]
z = vhel[galname]/3e5 
H0 = 70 #km/s/Mpc
d_prop = v/H0 #in Mpc
d_comov = d_prop * (1+z)
d = d_comov
pixelscale = 0.29*3 #arcsec/pix
#R = d*(0.29*3/2)/(360*60*60)*1e3 #in kpc
R = 2*math.pi*d*(pixelscale/3600)/(360)*10**6 #in pc

#Read in masking regions from GIST config files
mask = np.genfromtxt('../../configFiles/spectralMasking_PPXF', dtype = None,
                    skip_header = 3, names = ['lam','width','line'])
em_ndx = (mask['line'] != b'sky') #index for masked regions that are not sky lines

#Plot spectra from all bins. If you only want central bin, define nbins = [0]
fig, axes = plt.subplots(len(nbins),1)
fig.suptitle(folder)
for nbin in nbins:
    if len(nbins) > 1: 
        ax = axes.ravel()[nbin]
    else:
        ax = axes
    ax.plot(np.exp(lam), spec[nbin], 'k') #input spectra
    try: 
        ax.plot( np.exp(lam), fit['BESTFIT'][nbin], 'r', lw = 3) #best fit of continuum + emission lines
    except NameError:
        print("")
    ax.plot( np.exp(lam), cfit['BESTFIT'][nbin], 'orange', lw = 2) #continuum best fit
    #plt.text(4500,200+nbin*300, 'Bin '+str(nbin))
    if nbin == 0:
        ax.text(4500,1200, 'Bin 0 - Continuum center', fontsize = 20)
    else:
        ax.text(4500,1200, 'Bin '+str(nbin),fontsize = 20)#+' - 0.87" ('+str(int(R))+' pc) off center', 
        ax.text(4500,900, 'Sigma = '+str(sig_gist)[0:4]+"+/-"+str(sig_gist_err)[0:4],fontsize = 10) 
    ax.set_ylabel('Flux (in counts)', fontsize = 22)
    ax.set_xlim(np.exp(lam)[0],np.exp(lam)[-1])

    #plt.text(4500,775, 'Bin 2 - 0.87" ('+str(int(R))+' pc) off center', 
#         fontsize = 20)
    for mlam,width in zip(mask['lam'][em_ndx],mask['width'][em_ndx]):
        ax.axvspan(mlam-width/2,mlam+width/2,facecolor = 'k', alpha = 0.1)
    for mlam,width in zip(mask['lam'][~em_ndx],mask['width'][~em_ndx]):
        ax.axvspan(mlam-width/2,mlam+width/2,facecolor = 'blue', alpha = 0.1)
plt.xlabel('Wavelength (in Angstroms)', fontsize = 22)


###############################################################################
#Plot NII, SII and OI plots    
#Table with emission line fluxes
em = Table.read(folder+'_gas_BIN.fits')
#em = em[np.median(em, axis = 1) != -1.0]
#em_lam = Table.read(folder+'_gas_BIN.fits', hdu = 1)['_lambda']
em_name = em.keys()

#old naming convension - obsolete
#nii = em[:,em_name == '[NII]'][[0,2]]
#sii = em[:,em_lam == 6716.31][[0,2]]+ em[:,em_lam == 6730.68][[0,2]]#em[:,em_name == '[SII]']
#oi = em[:,em_name == '[OI]'][[0,2]]
#oiii = em[:,em_lam == 5006.77][[0,2]]#[:,em_name == '[OIII]']
#halpha = em[:,em_name == 'Ha'][[0,2]]
#hbeta = em[:,em_name == 'Hb'][[0,2]]

nii = em['[NII]_6583.34_F']
sii = em['[SII]_6716.31_F']+ em['[SII]_6730.68_F']#em[:,em_name == '[SII]']
oi = em['[OI]_6300.2_F']
oiii = em['[OIII]_5006.77_F']#[:,em_name == '[OIII]']
halpha = em['Ha_6562.8_F']
hbeta = em['Hb_4861.32_F']
#hbeta[0] = 2*hbeta[0]
error = 0
cmap = plt.get_cmap('rainbow', len(nii)*2)
#tag = np.arange(len(nii))
tag = np.array([0,1,1,3.4,8,3,7])*pixelscale #only for rs0010
n2ha = np.log10(nii/halpha)
s2ha = np.log10(sii/halpha)
o1ha= np.log10(oi/halpha)
o3hb = np.log10(oiii/hbeta)
tag = np.arange(len(n2ha))
if error: 
    n2ha_err = ratioerror(nii, nii_err, halpha, halpha_err)
#    s2ha_err = ratioerror(sii, sii_err, halpha, halpha_err)
    o1ha_err = ratioerror(oi, oi_err, halpha, halpha_err)
    o3hb_err = ratioerror(oiii, oiii_err, hbeta, hbeta_err)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.3,100)
cen = 0
nii_cen = nii[cen]
oiii_cen = oiii[cen]
#sii_cen = sii[cen]
oi_cen = oi[cen]
ha_cen = halpha[cen]
hb_cen = hbeta[cen]
fig, ax1 = plt.subplots()
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
cax = ax1.scatter(np.array(n2ha).ravel(), np.array(o3hb).ravel(), c= tag, cmap = cmap)
fig.colorbar(cax, extend='min')
#, alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.plot(np.log10(nii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', 
         markersize = 10,mfc = 'none', mew = 2)#, label = 'Definite Star Forming')
#ax1.set_xlim(-1.5,0.5)
#ax1.set_ylim(-1.0,1.0)
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax1.colorbar(extend='min')

if error: 
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')

#SII/OIII plot
#fig, ax2 = plt.subplots()
#ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
#ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
#                  'k--', label = 'Liner/Seyfert Division')
#cax = ax2.scatter(np.array(s2ha).ravel(), np.array(o3hb).ravel(), c= tag, cmap = cmap)
#fig.colorbar(cax, extend='min')
##, alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
#ax2.plot(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', 
#         markersize = 10,mfc = 'none', mew = 2)#, label = 'Definite Star Forming')
##ax2.set_xlim(-1.7, 1.7)
##ax2.set_ylim(-1.0,1.5)
#ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#if error:
#    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
#                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
#                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
#
#OI/OIII plot
o1ha[np.where(o1ha == -np.inf)] = -100
fig, ax3 = plt.subplots()
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                               'k--', label = 'Ke06 Liner/Seyfert Division Line')
ax3.scatter(np.array(o1ha).ravel(), np.array(o3hb).ravel(), c= tag, cmap = cmap)
fig.colorbar(cax, extend='min')
ax3.plot(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', 
         markersize = 10,mfc = 'none', mew = 2)#, label = 'Definite Star Forming')
ax3.set_xlim(-2.0, 0.5)
#ax3.set_ylim(-1.0,2.0)
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
if error:
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                                yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                            markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
