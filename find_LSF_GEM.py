import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import os
import glob
import pdb
'''
This code should extract relevant info in order to find the LSF from GEM data
files with sky spectra have the prefit 'steqxbprgS*' on them
multi extension fits, just load in the file and do fits.info()
the one label sky is the sky spectra 

NOTE: for blue setup data, we chose to drop the first chip worth of data for any processing because the FWHM is systematically larger. I am making sure that is also included in this as you don't want to fit that
'''

root = '/srv/two/sheila/derrcarr/GEM_DC/2013B/'
galname = 'rf0078'

skyfiles = glob.glob(root + galname + '/' + 'steqxbprg*') #prefix for post-sky subtraction file usually 2 sciences with slight 5 angstrom spec dither
arcfiles = glob.glob(root + galname + '/' + 'teprg*') #prefix for arc files


skylines = [5577.334, 6300.304, 6363.780, 6863.955] # red setup
arclines = [4702.3161,4806.0205,5062.0371,5090.4951] # blue setup

def extract_spectra_parameters(skyfiles,arcfiles,setup):
    if setup == 'red':        
        ### defining useful variables ###
        fitsfile = fits.open(skyfiles[0]) #opens first fits arbitrarily
        hdu = fitsfile[0].header
        skyspec = fitsfile[3].data # this should be the spectra of JUST the sky averaged over all sky fibers
        skyhdr = fitsfile[3].header # for some reason there are two sky fits for red setup, if this one is weird try [4]
        #pdb.set_trace()
        len_of_wavelength_axis = len(skyspec)
        delta_lamb = skyhdr['CDELT1']
        starting_lam = skyhdr['CRVAL1']
        lam = np.arange(starting_lam, starting_lam + delta_lamb * (len_of_wavelength_axis-1), delta_lamb)
        pix_to_fit = 30
        waverange =  [np.min(lam), np.max(lam)]
        conv=((waverange[1]-waverange[0])/(len_of_wavelength_axis-1))
        return lam, skyspec, conv, pix_to_fit
    if setup == 'blue': #for future reference, only applies to data before 2015. 
        ### defining useful variables ###
        fitsfile = fits.open(arcfiles[0]) #opens first arc fits arbitrarily
        hdu = fitsfile[0].header
        arcspec = fitsfile[2].data[0] # this should be the spectra of one random fiber
        archdr = fitsfile[2].header
        pdb.set_trace()
        len_of_wavelength_axis = len(arcspec)
        delta_lamb = archdr['CDELT1']
        starting_lam = archdr['CRVAL1']
        lam = np.arange(starting_lam, starting_lam + delta_lamb * (len_of_wavelength_axis-1), delta_lamb)
        pix_to_fit = 30

        #Throwing away first chip of data because it has systematically lower resolution
        new_sel = np.where(lam > 4540)
        lam = lam[new_sel]
        arcspec = arcspec[new_sel]
        len_of_wavelength_axis = len(arcspec)
        starting_lam = lam[0]
        
        waverange =  [np.min(lam), np.max(lam)]
        conv=((waverange[1]-waverange[0])/(len_of_wavelength_axis-1))
        return lam, arcspec, conv, pix_to_fit

#be sure to change name at the top
lam, skyspec, conv, pix_to_fit = extract_spectra_parameters(skyfiles,arcfiles,'red')

#plt.plot(lam,skyspec)
#plt.show()
