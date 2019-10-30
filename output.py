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

root = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/results/'
#galname = 'rf0376'
galname = 'rs0010' #strong OI
#galname = 'rf0503' #strong OI - very compact; only bins 0 and 2 have good SNR

#galname = 'rf0477' #- strong OI; good SNR ~ 7
#galname = 'rs0105' - strong OI but weak SNR for spectrum
#galname = 'rs1143'
#galname = 'rs0124' #very low SNR throughout : GIST coulg not bin upto SNR ~ 5
inputname = 'binned3d'+galname+'crop'
#inputname = 'binned3arcsdssrs0010'
#inputname = 'sami2sdssrs0010blue' #'sdss3arcrs0010'
#folder = root+inputname+'_1/'
folder = root+inputname+'_snr5/'
os.chdir(folder)

cfit = Table.read(inputname+'_ppxf-bestfit.fits')
fit = Table.read(inputname+'_gandalf-bestfit_BIN.fits')
data = fits.open('/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/'+inputname+'.fits')[0].data
hgal = fits.open('/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/'+inputname+'.fits')[0].header
lamdata = (np.arange(hgal['NAXIS3']) + hgal['CRPIX3']-1) * hgal['CDELT3'] + hgal['CRVAL3'] 
table = fits.open(folder+inputname+'_table.fits')[1].data
lam = fits.open(inputname+'_VorSpectra.fits')[2].data.LOGLAM
spec = fits.open(inputname+'_VorSpectra.fits')[1].data.SPEC

plam = fits.open(folder+inputname+'_ppxf-optimalTemplates.fits')[2].data['LOGLAM_TEMPLATE']
good = (lamdata > 4300) & (lamdata < 6900)
temp = fits.open(folder+inputname+'_ppxf-optimalTemplates.fits')
#plt.plot(np.exp(plam), np.transpose(temp[1].data[0]))
nbins = [0, 1, 2]
df= pd.read_csv('~/github/SDSS_spectra/RESOLVE_filter_new.csv')
df.index = df.name    

#db = readsav('/srv/one/resolve/database_internal/merged_idl_catalog/stable/resolvecatalog.dat')
 
db = readsav('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/resolvecatalog.dat')
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
mask = np.genfromtxt(folder+'spectralMasking_PPXF.config', dtype = None,
                    skip_header = 3, names = ['lam','width','line'])

plt.figure()
for nbin in nbins:
#    plt.plot(np.exp(lam), spec[nbin]+(nbin*2*np.max(spec[nbin-1])), 'k')
#    plt.plot( np.exp(lam), fit['BESTFIT'][nbin]+(nbin*2*np.max(spec[nbin-1])), 'r', lw = 3)
#    plt.plot( np.exp(lam), cfit['BESTFIT'][nbin]+(nbin*2*np.max(spec[nbin-1])), 'orange', lw = 2)
    plt.plot(np.exp(lam), spec[nbin]+(nbin*300), 'k')
#    plt.plot(lamdata/(1+z),data[:,:,23]*(1+z))
    plt.plot( np.exp(lam), fit['BESTFIT'][nbin]+(nbin*300), 'r', lw = 3)
    plt.plot( np.exp(lam), cfit['BESTFIT'][nbin]+(nbin*300), 'orange', lw = 2)
    #plt.text(4500,200+nbin*300, 'Bin '+str(nbin))

plt.text(4500,175, 'Bin 0 - Continuum center', fontsize = 20)
plt.text(4500,450, 'Bin 1 - 0.87" ('+str(int(R))+' pc) off center', 
         fontsize = 20)
plt.text(4500,775, 'Bin 2 - 0.87" ('+str(int(R))+' pc) off center', 
         fontsize = 20)
em_ndx = (mask['line'] != b'sky')
for mlam,width in zip(mask['lam'][em_ndx],mask['width'][em_ndx]):
    plt.axvspan(mlam-width/2,mlam+width/2,facecolor = 'k', alpha = 0.1)
for mlam,width in zip(mask['lam'][~em_ndx],mask['width'][~em_ndx]):
    plt.axvspan(mlam-width/2,mlam+width/2,facecolor = 'blue', alpha = 0.1)
plt.xlim(np.exp(lam)[0],np.exp(lam)[-1])
plt.xlabel('Wavelength (in Angstroms)', fontsize = 22)
plt.ylabel('Flux (in counts)', fontsize = 22)
#plt.plot(np.exp(plam), fit[0]['BESTFIT'])
#plt.plot(lamdata[good]/(1+z),data[:,:,24][good])#/max(data[:,:,24][good]))
#plt.plot(data[:,:,23][good])
#plt.plot(data[:,:,22][good])
    
em = Table.read(inputname+'_gandalf_BIN.fits', hdu = 2)['FLUX']
em = em[np.median(em, axis = 1) != -1.0]
em_lam = Table.read(inputname+'_gandalf_BIN.fits', hdu = 1)['_lambda']
em_name = Table.read(inputname+'_gandalf_BIN.fits', hdu = 1)['name']

#nii = em[:,em_name == '[NII]'][[0,2]]
#sii = em[:,em_lam == 6716.31][[0,2]]+ em[:,em_lam == 6730.68][[0,2]]#em[:,em_name == '[SII]']
#oi = em[:,em_name == '[OI]'][[0,2]]
#oiii = em[:,em_lam == 5006.77][[0,2]]#[:,em_name == '[OIII]']
#halpha = em[:,em_name == 'Ha'][[0,2]]
#hbeta = em[:,em_name == 'Hb'][[0,2]]
nii = em[:,em_name == '[NII]']
#sii = em[:,em_lam == 6716.31]+ em[:,em_lam == 6730.68]#em[:,em_name == '[SII]']
oi = em[:,em_name == '[OI]']
oiii = em[:,em_lam == 5006.77]#[:,em_name == '[OIII]']
halpha = em[:,em_name == 'Ha']
hbeta = em[:,em_name == 'Hb']
hbeta[0] = 2*hbeta[0]
error = 0
cmap = plt.get_cmap('rainbow', len(nii)*2)
#tag = np.arange(len(nii))
tag = np.array([0,1,1,3.4,8,3,7])*pixelscale #only for rs0010
n2ha = np.log10(nii/halpha)
#s2ha = np.log10(sii/halpha)
o1ha= np.log10(oi/halpha)
o3hb = np.log10(oiii/hbeta)
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
