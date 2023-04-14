#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 16:47:07 2020

@author: mugpol
This code should extract relevant info in order to find the LSF from SOAR and 
Gemini data. 
SOAR - files with sky spectra (i.e., without sky subtraction) have the prefix
       'lin_'. 'lincen_' files with the center slit are usually used.  
Gemini - files with sky spectra have the prefix 'steqxbprgS*'. They are multi 
        extension fits, just load in the file and do fits.info() the one label 
        sky is the sky spectra. 
        NOTE: for blue setup data, we chose to drop the first chip worth of 
        data for any processing because the FWHM is systematically larger. 
"""

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from mpfit import mpfit
import scipy.interpolate as interpolate
import os
import glob
import pandas as pd

def extract_spectra_soar(filename):
    '''
    Extract a 1-D spectrum from a given 2-D SOAR spectrum. Sky and arc files 
    are treated in the same manner. 
    '''
    
    specfile=fits.open(filename)
    spec=specfile[0].data
    hdu = specfile[0].header
    shap=np.shape(spec)
    
    #use a row with sky spectra that is off the galaxy - 100 rows above the 
    #center is arbitrarily chosen
    row = shap[0]//2+100
    spec = spec[row,:]

    #Wavelength solution of the galaxy
    '''
    NAXIS = Number of pixels on the x-axis (wavelength axis)
    CRPIX = Starting pixel number
    CDELT1 = delta(lambda) in Angstroms per pixel
    CRVAL1 = Wavelength value corresponding to starting pixel (CRPIX)
    '''    
    lam = ((np.arange(hdu['NAXIS1']) + hdu['CRPIX1']-1) * hdu['CDELT1']) + hdu['CRVAL1'] 

    return spec, lam

def lsf(spec,lam,lines,pix_to_fit):
    '''
    Find the Line Spread Function of an instrument by fitting gaussians to 
    the specified sky or arc lines using a 1-D sky/arc spectrum.
    '''
    
    waverange = [np.min(lam), np.max(lam)]  
    shap = np.shape(spec)
    conv=((waverange[1]-waverange[0])/(shap[0]-1)) #convolution kernel width

    lines=np.array(lines)
    pixlines=(lines/conv)-(waverange[0]/conv)
    pixlines = [int(x) for x in pixlines]
    x = []
#    x.append(np.arange(0,shap[0],1))
    for i in range(0,len(lines)):
        x.append(np.arange(pixlines[i]-pix_to_fit[i],pixlines[i]+pix_to_fit[i]))
    plt.figure()
    ax1=plt.subplot2grid((2,5),(0,0),colspan=5)
    
    plt.plot(lam,spec,'b-') 
    for i in range(0,len(lines)):
        plt.plot(lam[x[i]],spec[pixlines[i]-pix_to_fit[i]:pixlines[i]+pix_to_fit[i]],'r-')
    plt.title('red = arc lines used for FWHM calc')
    fwhminpix, fwhminA = measure_fwhm(spec,lines,pixlines,pix_to_fit,conv)
    return fwhminpix, fwhminA

def extract_spectra_gemini(specfiles,setup):
    '''
    Extract a 1-D spectrum from a given 2-D Gemini spectrum. Sky and arc files 
    are treated in the same manner. For blue setup data, we chose to drop the 
    data in the first chip because the FWHM is systematically larger. 
    There are no clear sky lines in blue setup spectra, so we use only arcs. 
    '''    
    if setup == 'red':        
        ### defining useful variables ###
        skyfiles = specfiles
        fitsfile = fits.open(skyfiles[0]) #opens first fits arbitrarily
        hdu = fitsfile[0].header
        skyspec = fitsfile[3].data # this should be the spectra of JUST the 
                                    #sky averaged over all sky fibers
        skyhdr = fitsfile[3].header # for some reason there are two sky fits 
                                #for red setup, if this one is weird try [4]
        
        len_of_wavelength_axis = len(skyspec)
        delta_lamb = skyhdr['CDELT1']
        starting_lam = skyhdr['CRVAL1']
        lam = np.arange(starting_lam, starting_lam + \
                        delta_lamb * (len_of_wavelength_axis-1), delta_lamb)
        return skyspec, lam

    if setup == 'blue': #for future reference, only applies to data before 2015. 
        ### defining useful variables ###
        arcfiles = specfiles
        fitsfile = fits.open(arcfiles[0]) #opens first arc fits arbitrarily
        hdu = fitsfile[0].header
        arcspec = fitsfile[2].data[0] # this should be the spectra of one random fiber
        archdr = fitsfile[2].header

        len_of_wavelength_axis = len(arcspec)
        delta_lamb = archdr['CDELT1']
        starting_lam = archdr['CRVAL1']
        lam = np.arange(starting_lam, starting_lam + \
                        delta_lamb * (len_of_wavelength_axis-1), delta_lamb)

        #Throwing away first chip of data because it has systematically lower resolution
        new_sel = np.where(lam > 4540)
        lam = lam[new_sel]
        arcspec = arcspec[new_sel]
        len_of_wavelength_axis = len(arcspec)
        starting_lam = lam[0]
        
        return arcspec, lam

def myfunct(p, fjac=None, x=None, y=None, err=None):
        model = p[0] * np.exp(-((x-p[1])**2.)/(2.*p[2]**2.)) + p[3]
        status = 0
        return([status, (y-model)/err])
def func2(x, a, b, d, c):
        return a * np.exp(-((x-b)**2.)/(2.*d**2.)) + c
 
def measure_fwhm(spec,lines,pixlines,pix_to_fit, conv):
    fwhminpix=[]
    for i in np.arange(len(lines)):
        nput=spec[pixlines[i]-pix_to_fit[i]:pixlines[i]+pix_to_fit[i]]
        gerr=np.zeros(len(nput))+0.5
        xaxis=np.arange(len(nput))
    
        max=np.max(nput)
        mid=np.argmax(nput)
        #print(i,mid,max)
    
        p0=[max,mid,5.,30.]

    
        fa = {'x':xaxis, 'y':nput, 'err':gerr}
    
        m=mpfit(myfunct,p0,functkw=fa)
    
        fitdata=func2(xaxis,m.params[0],m.params[1],m.params[2],m.params[3])
    
        #figg1=ax7.pix_to_fit_subplot(3,1,i+1)
        sigma=m.params[2]
        fwhm=np.abs(sigma*2.35482)
        fwhminpix.append(fwhm)
    
        ax2=plt.subplot2grid((2,5),(1,i))#,colspan=4)
        plt.plot(xaxis,nput)
        plt.plot(xaxis,fitdata)
        plt.title('gaussian fits to arc lines')
        plt.ylim([np.min(nput),np.max(nput)])
        
        ax2.annotate('FWHM =',xy=(np.min(xaxis)+3.,np.max(nput)/2.0))
        ax2.annotate('%3.3g' % (fwhm*conv) + ' A',xy=(np.min(xaxis)+4.,np.max(nput)/2.3))
    fwhminpix = np.array(fwhminpix)
    fwhminA=fwhminpix*conv
        
    return fwhminpix, fwhminA

###############################################################################
tel = 'soar'
tel = 'gemini'
if tel == 'soar':
    #Specifying input directories
    root = '/srv/two/resolve/working/spec/'
    galname = 'rf0477'
    #galname = 'rs0010'
    #galname = 'rf0376'
    folderlist = {'rs0010': 'DS_2012-04-15/',
              'rs0376': 'EH_2014-09-20/',
              'rf0477': 'KE_2012-07-26/'}
    folder = root+folderlist[galname]
    
    #Specify skylines in SOAR broadspec range
    skylines = [5577.334, 6300.304, 6363.780, 6863.955] 
    #Specify comp lines for SOAR broadspec
    #lines=[3650,4358,5460,6965] 
    comp_file = folder+'lincen_'+galname+'bspec.fits' #using gal spectra before 
                                                #sky subtraction to use skylines
    #comp_file = folder+"cen_tbz2718.HgAr.fits" #using arcfile 

    #number of pixels around the lines that should be fit
    pix_to_fit=[30, 30, 25, 10] 

    #Extracting a 1-D sky/comp SOAR spectrum 
    spec,lam = extract_spectra_soar(comp_file)

    #Computing the Line Spread Function for SOAR
    fwhminpix, fwhminA = lsf(spec,lam,skylines,pix_to_fit)    

    #Saving the output in a GIST-friendly format
    save = 0
    outputfile = 'lsf_SOAR'
    if save:
        np.savetxt(outputfile,list(zip(skylines,fwhminA)))

if tel == 'gemini':
    #Specifying input directories
    root = ''
    galname = 'rf0078'
    
    #Specify skylines in Gemini blue setup range
    arclines = [4702.3161,4806.0205,5090.4951] # blue setup 
    #5062.0371,
    #allarclines = pd.read_csv("CuAr_GMOS.dat", delimiter = "  ", comment = '#', \
    #                       names = ['lambda','comments']) #file with all arclines in Gemini blue setup
    #arclines = np.array(allarclines['lambda'])
    #arclines = arclines[(arclines >= waverange[0]) & (arclines <= waverange[1])]
    arcfiles = glob.glob(root + galname + '/' + 'teprg*') #prefix for arc files

    #Specify skylines in Gemini red setup range
    skylines = [5577.334, 6300.304, 6363.780, 6863.955] # red setup
    skyfiles = glob.glob(root + galname + '/' + 'steqxbprg*') #prefix for post-sky subtraction file usually 2 sciences with slight 5 angstrom spec dither
        
    specfile = arcfiles
    spec, lam = extract_spectra_gemini(specfile,'blue')
    
    #number of pixels around the lines that should be fit
    pix_to_fit=[15, 20, 5] 

    fwhminpix, fwhminA = lsf(spec,lam,arclines,pix_to_fit)    

    #Saving the output in a GIST-friendly format
    save = 0
    outputfile = 'lsf_GEMINI'
    if save:
        np.savetxt(outputfile,list(zip(skylines,fwhminA)))

