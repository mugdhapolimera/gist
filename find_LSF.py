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

comp_file = folder+'lincen_'+galname+'bspec.fits' #using gal spectra before sky subtraction to use skylines
arc=fits.open(comp_file)
arcspec=arc[0].data
hdu = arc[0].header
shap=np.shape(arcspec)

lines = [5577.334, 6300.304, 6363.780, 6863.955] #skylines in SOAR broadspec range
#lines=[3650,4358,5460,6965] #comp lines for SOAR broadspec
lam = ((np.arange(hdu['NAXIS1']) + hdu['CRPIX1']-1) * hdu['CDELT1']) + hdu['CRVAL1'] 
waverange = [np.min(lam), np.max(lam)]
add=[30, 30, 25, 10] #number of pixels around the line that are included in the fit
#add=[20, 20, 20, 20]


conv=((waverange[1]-waverange[0])/(shap[1]-1))

lines=np.array(lines)
pixlines=(lines/conv)-(waverange[0]/conv)
pixlines = [int(x) for x in pixlines]

x0=np.arange(0,shap[1],1)
x1=np.arange(pixlines[0]-add[0],pixlines[0]+add[0],1)
x2=np.arange(pixlines[1]-add[1],pixlines[1]+add[1],1)
x3=np.arange(pixlines[2]-add[2],pixlines[2]+add[2],1)
x4=np.arange(pixlines[3]-add[3],pixlines[3]+add[3],1)
#x5=np.arange(pixlines[4]-add[4],pixlines[4]+add[4],1)

plt.figure()
ax4=plt.subplot2grid((2,5),(0,0),colspan=5)
#use a row that is off the galaxy - 100 rows above the center is arbitrarily chosen
row = shap[0]//2+100
plt.plot(lam,arcspec[row,:],'b-') 

plt.plot(lam[x1],arcspec[row,pixlines[0]-add[0]:pixlines[0]+add[0]],'r-')
plt.plot(lam[x2],arcspec[row,pixlines[1]-add[1]:pixlines[1]+add[1]],'r-')
plt.plot(lam[x3],arcspec[row,pixlines[2]-add[2]:pixlines[2]+add[2]],'r-')
plt.plot(lam[x4],arcspec[row,pixlines[3]-add[3]:pixlines[3]+add[3]],'r-')
#plt.plot(lam[x5],arcspec[shap[0]//2+100,pixlines[4]-add[4]:pixlines[4]+add[4]],'r-')


plt.title('red = arc lines used for FWHM calc')

fwhminpix=[]
for i in np.arange(len(lines)):
    nput=arcspec[row,pixlines[i]-add[i]:pixlines[i]+add[i]]
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

    #figg1=ax7.add_subplot(3,1,i+1)
    sigma=m.params[2]
    fwhm=np.abs(sigma*2.35482)
    fwhminpix.append(fwhm)
    
    ax5=plt.subplot2grid((2,5),(1,i))#,colspan=4)
    plt.plot(xaxis,nput)
    plt.plot(xaxis,fitdata)
    plt.title('gaussian fits to arc lines')
    plt.ylim([np.min(nput),np.max(nput)])
    
    ax5.annotate('FWHM =',xy=(np.min(xaxis)+3.,np.max(nput)/2.0))
    ax5.annotate('%3.3g' % (fwhm*conv) + ' A',xy=(np.min(xaxis)+4.,np.max(nput)/2.3))

fwhminpix=np.array(fwhminpix)

fwhminA=fwhminpix*conv

np.savetxt('lsf_SOAR',list(zip(lines,fwhminA)))
