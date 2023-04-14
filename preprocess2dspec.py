'''
July 15th 2019

Preprocess SOAR spectral data before using as input to GIST pipeline. 
Step 1. Bin SOAR spectra by 3 rows starting in the center
Step 2. Cast 2D spectra into 3D cube 

WARNING: This code works is written using Python 3.6. It may not work as 
        intended on Python

@author: mugdhapolimera
'''

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from mpfit import mpfit
import scipy.interpolate as interpolate
import os

def findcenter(gal_lin):
    '''
    Function to find the center of continuum of the galaxy 
    by fitting a Gaussian to the continuum
    '''
    
    shape=np.shape(gal_lin)
    center=int(shape[0]/2.)
    
    box=[center-150,center+150,800,1000]
    sumforgauss=np.sum(gal_lin[box[0]:box[1],box[2]:box[3]],axis=1)
    offset=center-150.
    xaxis=np.arange(len(sumforgauss))+offset
    gerr=np.zeros(len(sumforgauss))+0.5
    
    amp=np.max(sumforgauss)
    cen=np.argmax(sumforgauss)+offset
    
    #run mpfit
    p0=[amp,cen,10.,30.]
    
    def myfunct(p, fjac=None, x=None, y=None, err=None):
        model = p[0] * np.exp(-((x-p[1])**2.)/(2.*p[2]**2.)) + p[3]
        status = 0
        return([status, (y-model)/err])
    
    fa = {'x':xaxis, 'y':sumforgauss, 'err':gerr}
    
    m=mpfit(myfunct,p0,functkw=fa)
    #print (m.params[1])
    return int(m.params[1])

def bin_spectra(filename, bin_type, crop, standardflag,
                standard_name = '/srv/one/resolve/analysis_code/emlineanalysis/broademline/R2013-11-04_LTT1788_cen_broad_sens.fits'):
    """
    filename (str) = Filename of 2-D galaxy spectra
    compfilename (str) = Filename of 2-D arc file
    standardflag (0/1) = Use Standard Star spectrum for relative flux calibration    
    bin_type (str) = 'rows' for binning 3 rows starting from center
                     'sdss' for mimicing SDSS 3" aperture 
    standard_name (str) = Filename of standard star spectrum    
    """
    galaxy = fits.open(filename)
    gal = galaxy[0].data
    hgal = galaxy[0].header

    #Wavelength solution of the galaxy
    #NAXIS = Number of pixels on the x-axis (wavelength axis)
    #CRPIX = Starting pixel number
    #CDELT1 = delta(lambda) in Angstroms per pixel
    #CRVAL1 = Wavelength value corresponding to starting pixel (CRPIX)
    lam = (np.arange(hgal['NAXIS1']) + hgal['CRPIX1']-1) * hgal['CDELT1'] + hgal['CRVAL1'] 

    if standardflag:
        standard = fits.open(standard_name)
        std = standard[0].data
        hstd = standard[0].header
    
        #Wavelength solution of the standard star
        stdlam = (np.arange(hstd['NAXIS1']) + (hstd['CRPIX1']-1))*hstd['CDELT1'] + hstd['CRVAL1'] 
        
        #Interpolate flux calibration to lam (galaxy spectrum wavelengths)
        #The wavelength solution of the Standard Star and Galaxy may be differernt, 
        #so resample the standard star flux to match the galaxy's wavelengths
        f = interpolate.interp1d(stdlam,std, fill_value = 'extrapolate')
        flux = f(lam)
        fluxcorr = 10**(-0.4*(flux - np.median(flux))) #if sens in mag, gives ratios
    
        #Multiply galaxy spectrum by fluxcorr for relative calibration (not absolute)
        galcorr = gal*fluxcorr

        #Plot to check
        #plt.figure()
        #plt.plot(stdlam, std,'b')
        #plt.plot(lam,flux,'r')
        #plt.plot(lam,gal[int(len(gal)/2)],'g')
        #plt.plot(lam,galcorr[int(len(gal)/2)],'k')
        gal = galcorr
    
    
    #Find the center of the galaxy
    #Quick and dirty way to find the center
    #center = np.where(gal[:, int(np.shape(gal)[1]/2)] == \
    #                   max(gal[:, int(np.shape(gal)[1]/2)]))[0]
    
    #Finding center of the galaxy by finding mean of gaussian fit to the continuum
    center = findcenter(galcorr)

    if bin_type == 'rows':
        #If the number of rows in the galaxy is not a multiple of 6, zero pad rows to
        #the closest multiple of 6.
        #This is done so that the center of the new binned array is an integer number
        shape = np.shape(galcorr)
        if shape[0]%6:
            rows = 6 - (shape[0]%6)
            galcorr = np.append(galcorr,np.zeros([rows,shape[1]]), axis = 0)
        shape = np.shape(galcorr)
        
        #Create an empty array to store the binned galaxy with 1/3rd the number of rows
        binned_gal = np.zeros([shape[0]//3, shape[1]])
        cen_bin = np.shape(binned_gal)[0]//2 #center of the new binned spectrum

        '''
        Binning from the center of the galaxy; every row in the binned galaxy gets the
        value of the total flux in 3 rows of the original flux array. This is done
        since there is no independent information in 3 adjacent rows of the original
        spectrum.
        '''
        for i in np.arange(0,np.shape(binned_gal)[0]//2):
            binned_gal[cen_bin+i] = np.sum(galcorr[(center+(3*i) -1):\
                      (center+(3*i) +1)], axis = 0)
        
            binned_gal[cen_bin-i] = np.sum(galcorr[(center-(3*i) -1):\
                      (center-(3*i) +1)], axis = 0)
        
        if crop:
        #Trim the new binned spectra to have fewer non-galaxy spaxels
            binned_gal = binned_gal[cen_bin-(binned_gal.shape[0]//5): \
                                cen_bin+(binned_gal.shape[0]//5)]
    if bin_type == 'sdss':
        #Binning 2 arcsecs around the center 
        binned_gal = np.zeros([1, shape[1]])
        binned_gal[0,:] = np.sum(galcorr[(center-3):(center+3)],axis = 0)
    #binned_gal = binned_gal[None,:]
    return binned_gal

def make_cube(spec2d, hdu, outputfile):
    #Add a dummy dimension and make the longslit spectrum a 3D data-cube
    #Add keywords to the header corresponding to the 3rd dimension as wavelength
    spec3d = spec2d.T[:,None]
    hdu['NAXIS'] = 3
    hdu['NAXIS1'] = spec3d.shape[1]
    hdu['NAXIS2'] = spec3d.shape[2]
    hdu['NAXIS3'] = spec3d.shape[0]
    
    hdu['CRVAL3'] = hdu['CRVAL1']
    hdu['CRVAL1'] = 1
    
    hdu['CRPIX3'] = 1
    
    hdu['CDELT3'] = hdu['CDELT1']
    hdu['CD3_3'] = hdu['CD1_1']
    hdu['CDELT1'] = 1
    hdu['CD1_1'] = 1
       
    if not os.path.exists(outputfile):
        fits.writeto(outputfile, spec3d, hdu)
    return spec3d


root = '/srv/two/resolve/working/spec/'
galname = 'rf0477'
folderlist = {'rs0010': 'DS_2012-04-15/',
          'rs0346': 'EH_2014-09-20/',
          'rf0477': 'KE_2012-07-26/'}
folder = root+folderlist[galname]
filename = folder+'slincen_'+galname+'bspec.fits'
comp_file = folder+'cen_arc_broad.fits'
bin_type = 'rows'; crop = 1; standardflag = 1
binnedspec2d = bin_spectra(filename, bin_type, crop, standardflag)
outputfile = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/newgisttutorial/gistTutorial/inputData/binned3d'+galname+'crop.fits'
#outputfile = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/sdss2arc'+galname+'.fits'
hdu = fits.getheader(filename)
spec3d = make_cube(binnedspec2d, hdu, outputfile)
