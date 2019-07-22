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





df = pd.read_csv('../SDSS_spectra/RESOLVE_bpt1filter_new.csv')
df.index = df.name
RunMC = True


galname = 'rf0376'
filename = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/binned3drf0376crop.fits'

h1 = fits.open(filename)[0].header
summedspec=fits.read(filename)[0].data[2:-5]

#Program reads off the wavelength range for this spectrum
wave_md = h1['CRVAL3']
wave_pix = h1['CRPIX3']
wave_del = h1['CD3_3']
wave_st = wave_md - wave_del*wave_pix + 2*wave_del
wave_ed = wave_md + wave_del*wave_pix - 5*wave_del
waves = [wave_st, wave_ed]

print('')
print('Starting Wavelength is ' + str(wave_st) + ' Angsroms')
print('Ending Wavelength is ' + str(wave_ed) + ' Angsroms')
print('')

lam = np.linspace(waves[0], waves[1], num=len(summedspec))
lamRange1=np.array([waves[0], waves[1]])	#start wavelength, end wavelength

## Truncate off low resolution part of our data
## Blue region to the left of the 1st chip gap FWHM = 1.2
## Rest of the image further to the red FWHM = 1.0
if wave_st < 4070:
    gap = 4540 #chipgap blue
else:
    gap = 4580 #chipgap red
sel = np.where(lam > gap)
summedspec = summedspec[sel]
lam = lam[sel]
waves[0] = np.min(lam)
lamRange1[0] = np.min(lam)



lamRange1, lam, = read_gemini(specfile)

galaxy, logLam1, velscale = util.log_rebin(lamRange1, summedspec)
median = np.nanmedian(galaxy)
galaxy = galaxy / median
noise = np.ones_like(galaxy) * h1['NOISE'] / median #constant noise based on 
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
z2 = copy.copy(z)
z=0  #Bring Galaxy to Rest Frame


sspNew, logLam2, velscale_temp = util.log_rebin(lamRange2, templates[:,0],
                                        velscale=velscale/velscale_ratio)
templatesNew = np.zeros((sspNew.shape[0], templates.shape[1]))

for s in range(templates.shape[1]):
        ssp = templates[:,s]
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
    galaxy = galaxy + normal


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