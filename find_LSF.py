#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 16:47:07 2020

@author: mugpol
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 16:29:55 2020

@author: mugdhapolimera
"""

'''
July 15th 2019

Script to preprocess SOAR spectral data before using as input to gist pipeline

@author: mugdhapolimera
'''

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from mpfit import mpfit
import scipy.interpolate as interpolate
import os

root = '/srv/two/resolve/working/spec/'
galname = 'rf0477'
#galname = 'rs0010'
#galname = 'rf0376'
folderlist = {'rs0010': 'DS_2012-04-15/',
          'rs0376': 'EH_2014-09-20/',
          'rf0477': 'KE_2012-07-26/'}
folder = folderlist[galname]
folder = root+folder

#Wavelength solution of the galaxy
#NAXIS = Number of pixels on the x-axis (wavelength axis)
#CRPIX = Starting pixel number
#CDELT1 = delta(lambda) in Angstroms per pixel
#CRVAL1 = Wavelength value corresponding to starting pixel (CRPIX)

lines = [5577.334, 6300.304, 6363.780, 6863.955] #skylines in SOAR broadspec range
#lines=[3650,4358,5460,6965] #comp lines for SOAR broadspec
comp_file = folder+'lincen_'+galname+'bspec.fits' #using gal spectra before sky subtraction to use skylines
pix_to_fit=[30, 30, 25, 10] #number of pixels around the line that are included in the fit
fwhminpix, fwhminA = lsf(comp_file,lines,pix_to_fit)    
def lsf(filename,lines,pix_to_fit):
    arc=fits.open(filename)
    spec=arc[0].data
    hdu = arc[0].header
    shap=np.shape(spec)
    lam = ((np.arange(hdu['NAXIS1']) + hdu['CRPIX1']-1) * hdu['CDELT1']) + hdu['CRVAL1'] 
    waverange = [np.min(lam), np.max(lam)]  

    conv=((waverange[1]-waverange[0])/(shap[1]-1)) #convolution kernel width

    lines=np.array(lines)
    pixlines=(lines/conv)-(waverange[0]/conv)
    pixlines = [int(x) for x in pixlines]

    x = []
    x.append(np.arange(0,shap[1],1))
    for i in range(1,len(lines)):
        x.append(np.arange(pixlines[i]-pix_to_fit[i],pixlines[i]+pix_to_fit[i]))

    plt.figure()
    ax4=plt.subplot2grid((2,5),(0,0),colspan=5)
    
    #use a row that is off the galaxy - 100 rows above the center is arbitrarily chosen
    row = shap[0]//2+100
    spec = spec[row,:]
    plt.plot(lam,spec,'b-') 
    for i in range(1,len(lines)):
        plt.plot(lam[x[i]],spec[pixlines[i]-pix_to_fit[i]:pixlines[i]+pix_to_fit[i]],'r-')
    plt.title('red = arc lines used for FWHM calc')
    fwhminpix, fwhminA = measure_fwhm(spec,lines,pixlines,pix_to_fit)
    return fwhminpix, fwhminA
 
def measure_fwhm(spec,lines,pixlines,pix_to_fit):
    fwhminpix=[]
    for i in np.arange(len(lines)):
        nput=spec[pixlines[i]-pix_to_fit[i]:pixlines[i]+pix_to_fit[i]]
        gerr=np.zeros(len(nput))+0.5
        xaxis=np.arange(len(nput))
    
        max=np.max(nput)
        mid=np.argmax(nput)
        print(i,mid,max)
    
        p0=[max,mid,5.,30.]

    def myfunct(p, fjac=None, x=None, y=None, err=None):
        model = p[0] * np.exp(-((x-p[1])**2.)/(2.*p[2]**2.)) + p[3]
        status = 0
        return([status, (y-model)/err])

    fa = {'x':xaxis, 'y':nput, 'err':gerr}

    m=mpfit(myfunct,p0,functkw=fa)

    def func2(x, a, b, d, c):
        return a * np.exp(-((x-b)**2.)/(2.*d**2.)) + c

    fitdata=func2(xaxis,m.params[0],m.params[1],m.params[2],m.params[3])

    #figg1=ax7.pix_to_fit_subplot(3,1,i+1)
    sigma=m.params[2]
    fwhm=np.abs(sigma*2.35482)
    fwhminpix.append(fwhm)
    fwhminpix = np.array(fwhminpix)
    fwhminA=fwhminpix*conv
    
    ax5=plt.subplot2grid((2,5),(1,i))#,colspan=4)
    plt.plot(xaxis,nput)
    plt.plot(xaxis,fitdata)
    plt.title('gaussian fits to arc lines')
    plt.ylim([np.min(nput),np.max(nput)])
    
    ax5.annotate('FWHM =',xy=(np.min(xaxis)+3.,np.max(nput)/2.0))
    ax5.annotate('%3.3g' % (fwhm*conv) + ' A',xy=(np.min(xaxis)+4.,np.max(nput)/2.3))
    
    return fwhminpix, fwhminA


np.savetxt('lsf_SOAR',list(zip(lines,fwhminA)))
