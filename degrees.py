#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 11:46:49 2019

Testing which Additive and Multiplicative Legendre Polynomials are the 
optimal for pPXF fitting. 

Adapted from runppxf.py by Nate Cleaves

@author: mugpol
"""

import numpy as np
import time
import copy
import pandas as pd
from astropy.io import fits

import multiprocessing as mp
import matplotlib.pyplot as plt

import glob

from ppxf import ppxf
from ppxf import ppxf_util as util
from scipy import ndimage
model = 'maraston'

# =====================================================================================

def DER_SNR(flux):
   
# =====================================================================================
   from numpy import array, where, median, abs 

   flux = array(flux)

   # Values that are exactly zero (padded) are skipped
   flux = array(flux[where(flux != 0.0)])
   n    = len(flux)      

   # For spectra shorter than this, no value can be returned
   if (n>4):
      signal = median(flux)

      noise  = 0.6052697 * median(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
      return noise #float(signal / noise)  

   else:

      return 0.0

# end DER_SNR -------------------------------------------------------------------------

def read_models(vazdekis, z):


    test = fits.open(vazdekis[0])
    sh = test[0].data.shape
    test.close()

    #loads in models

    templates = np.zeros((sh[0], len(vazdekis)))
    for j in range(len(vazdekis)):
        if (j % 5 == 0):
            print(str(j) +  ' / ' + str(len(vazdekis)))
        hdu = fits.open(vazdekis[j])
        ssp = hdu[0].data
        h2 = hdu[0].header
        print (h2)
        if model != 'kurucz':
            lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
        else:
            lam2 = hdu[1].data
            lamRange2 = [np.min(lam2), np.max(lam2)]
        try:
            h2['CRVAL1']
        except:
            h2['CRVAL1'] = np.min(lam2)

        if model == 'kurucz':
            #Linear resampling of kurucz data
            #Assumes nearly linear sampling in the area of interest

            ssp, tdl = resample(ssp, lam2, np.mean(lam2[1:] - lam2[:-1]))
            h2['CDELT1'] = tdl
            lamRange2 = lam2[0] + np.array([0.,tdl*(h2['NAXIS1']-1)])
            lam2 = np.linspace(lamRange2[0], lamRange2[1], num=h2['NAXIS1'])

        templates[:,j] = ssp
        hdu.close()
    plt.plot(templates)

    return h2, lamRange2, templates

def minimize_leg(galname, templates, galaxy, noise, velscale, start, goodPixels,
                 dv, velscale_ratio, iterations, useMC, find_min=True):

#Minimizes ppxf uncertainty using legendre polynomials
#find_min=True will force the program to find the minimum again
    value_range = range(5,23)

    if useMC == False: #so that you don't double up on the multiprocessing and kill the server
        #iterations
        feeds = []
        for ad in value_range:
            for md in value_range:
                feeds.append([(ad, md), i, iterations, templates, galaxy, 
                noise, velscale, start, goodPixels, dv, velscale_ratio])
        P = mp.Pool(iterations)
        outputs = P.map(min_worker, feeds)

        outputs = np.array(outputs)
        degrees = outputs[:,0]
        errors = outputs[:,1]
        deg = degrees[np.argmin(errors)]
        

        print(' ')
        print('Deg is:')
        print(deg)
        print(' ')
        
        return deg
        

    
    else: #serial version
        
        stime = time.time()
        
        deg_used = []
        errors   = []

        if find_min or np.any(np.isnan(degrees)):
        #Sample parameter space within $value_range to find minimum uncertainty
            minimum = 1e10
            degrees = [0,0] #ad, md
            for ad in value_range: #degree of additive legendre polynomial
                for md in value_range: #degree of multiplicative legendre polynomial
                    print('degrees of %i additive, %i multiplicative'%(ad, md))
                    #print((ad, md))
                    
                    ppxf_obj = ppxf.ppxf(templates, galaxy, noise, velscale, 
                                         start, goodpixels=goodPixels, 
                                         plot=False, moments=4, degree=ad,
                                         vsyst=dv, 
                                         velscale_ratio = velscale_ratio, 
                                         mdegree=md, clean=False)
                    error = ppxf_obj.error[1]*np.sqrt(ppxf_obj.chi2)
                    
                    deg_used.append((ad,md))
                    errors.append(error)
                    
                    #error = ppxf_obj.chi2
                    print('Uncertainty is ' + str(round(float(error), 2)) + ' km/s')
                    if error < minimum:
                        minimum = error
                        degrees[0] = ad
                        degrees[1] = md
                
        print('')
        print('Degree of Linear Legendre Polynomial is: ' + str(degrees[0]))
        print('')
        print('Degree of Multiplicative Legendre Polynomial is: ' + str(degrees[1]))
        print('')
        print('Took %s seconds to complete'%int(round(time.time() - stime)))
        print('')


    return degrees



iterations = 30
models = ['maraston', 'vazdeki', 'miles', 'elodie', 'kurucz']

model = 'maraston' #ask_user_options('Which models should we use?', models)
if model == 'maraston':
        #SSPs made based on Prof. Maraston's work by Mark Norris
    vazdekis=glob.glob('/afs/cas.unc.edu/users/m/a/manorris/public/kroupa/templates/*.fits')
    #vazdekis=glob.glob('/srv/two/sheila/natetc/ppxfblue/marastontest/*.fits')
    FWHM_tem = 0.55
elif model == 'vazdeki':
        #Stock vazdeki stellar library provided with ppxf
    vazdekis = glob.glob(ppxf_dir + '/miles_models/Mun1.30Z*.fits')
    FWHM_tem = 2.51
elif model == 'elodie':
        #Elodie stellar library
    vazdekis = glob.glob('/afs/cas.unc.edu/users/n/a/natetc/public/elodie/H/*.fits')
    FWHM_tem = 0.5
elif model == 'miles':
        #Miles stellar library
    vazdekis = glob.glob('/afs/cas.unc.edu/users/n/a/natetc/public/miles/*.fits')
    FWHM_tem = 2.51
elif model == 'kurucz':
        #Model stellar spectra produced by Kurucz; basically a very small stellar library
    vazdekis = glob.glob('/srv/two/sheila/natetc/ppxfblue/kurucz*.fits')
    FWHM_tem = []
    for vazdeki in vazdekis:
        vaz = fits.open(vazdeki)
        FWHM_tem.append(vaz[0].header['FWHM'])

df = pd.read_csv('../SDSS_spectra/RESOLVE_bpt1filter_new.csv')
#df = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_spectra/RESOLVE_filter_new.csv')

df.index = df.name
RunMC = True


galname = 'rs0010'
filename = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/binned3drs0010crop.fits'
#filename = 'binned3drf0376crop.fits'

h1 = fits.open(filename)[0].header
spec=fits.open(filename)[0].data
summedspec = np.sum(spec, axis = 2)[:,0]

#Program reads off the wavelength range for this spectrum

lam = (np.arange(h1['NAXIS3']) + h1['CRPIX3']-1) * h1['CD3_3'] + h1['CRVAL3']
lamRange1=np.array([lam[0], lam[-1]])	#start wavelength, end wavelength

print('')
print('Starting Wavelength is ' + str(lam[0]) + ' Angstroms')
print('Ending Wavelength is ' + str(lam[-1]) + ' Angstroms')
print('')




#error calculated around continuum by binning program
FWHM_gal = 6.5 #from sky line fitting
c = 299792.458
vel = df.loc[galname].cz
z = vel/c
velscale_ratio = 1

h2, lamRange2, templates = read_models(vazdekis, z)


### Cosmological Redshift Correction ###    
lamRange1 = lamRange1 / (1+z)
FWHM_gal = FWHM_gal / (1+z)

goodlam = (lam/(1+z) > 3900) & (lam/(1+z) < 6797)
lam = lam[goodlam]
lamRange1 = np.array([lam[0], lam[-1]]) / (1+z)
summedspec = summedspec[goodlam]

gal, logLam1, velscale = util.log_rebin(lamRange1, summedspec)
median = np.nanmedian(gal)
#galaxy = galaxy / median
noise = np.zeros(spec.shape[0])
for i in range(spec.shape[0]):
    noise[i] = DER_SNR(spec[i,:,:])
noise = noise[goodlam]/median

z2 = copy.copy(z)
z=0  #Bring Galaxy to Rest Frame
gal = gal/median
#galaxy = galaxy[overlaplam]
#noise = noise[overlaplam]
lam_temp = h2['CRVAL1'] + h2['CDELT1']*np.arange(h2['NAXIS1'])
    

FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
sigma = FWHM_dif/2.355/h2['CDELT1'] #dl/pix
#print(sigma)
blur_temp = True

sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, templates[:,0],
                                        velscale=velscale/velscale_ratio)
templatesNew = np.zeros((sspNew.shape[0], templates.shape[1]))

for s in range(templates.shape[1]):
        ssp = templates[:,s]
        if blur_temp:
            ssp = ndimage.gaussian_filter1d(ssp,sigma)
        sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, ssp,velscale=velscale/velscale_ratio)
        templatesNew[:,s] = sspNew/np.median(sspNew)

    
templates = templatesNew

if velscale_ratio > 1:
    dv = (np.mean(logLam2[:velscale_ratio]) - logLam1[0])*c  # km/s
else:
    dv = (logLam2[0] - logLam1[0])*c  # km/s

sig = 3*velscale
    
start = [vel, sig]

stime = time.time()

if RunMC:
    microsecond = int(str(time.time()).split('.')[-1])
    #Randomize numpy seeds. Otherwise all processes will give exact same values. 
    #This seed depends on exact system time
    np.random.seed(seed=microsecond) 
    normal = np.random.normal(size=noise.shape) * noise
    galaxy = gal + normal

goodPixels = np.arange(len(galaxy))

em_lam = np.array([4363, 5007, 4861, 6548, 6562, 6584, 6717, 6731])

goodPixels = [x for x in range(len(lam)) if ~((lam[x] < (em_lam + 10)) & (lam[x] > (em_lam - 10))).all()]
goodPixels = np.array(goodPixels)
degrees = [10,10]
degrees = minimize_leg(galname, templates, galaxy, noise, velscale, start, 
                      goodPixels, dv, velscale_ratio, iterations, RunMC, 
                     find_min = True)

ppxf_obj = ppxf.ppxf(templates, galaxy, noise, velscale, start, 
                     goodpixels=goodPixels, plot=False, moments=4,
                     degree=degrees[0],vsyst=dv, 
                     velscale_ratio = velscale_ratio, 
                     mdegree=degrees[1], clean=False)

redshift = ppxf_obj.sol[0]
vel_disp = ppxf_obj.sol[1]
vel_disp_err = ppxf_obj.error[1]*np.sqrt(ppxf_obj.chi2)

print("Minimization took %s minutes to complete"
      %int(round((time.time() - stime) / 60.)))

ppxf_obj.plot()