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
root = '/srv/two/resolve/working/spec/'
#folder = root+'EH_2012-03-15/' #rs0105
#folder = 'DR_2012-03-19/' #rs1143
#folder = 'EH_2012-03-26/' #rs0124
#folder = 'KE_2012-05-18/' #rs1283
#folder = 'DS_2012-10-10/' #rf0503
#filename = folder+'slincen_rf0376bspec.fits'
galname = 'rf0477'
#galname = 'rs0010'
#galname = 'rs0105'
#galname = 'rs1143'
#galname = 'rs0124'
#galname = 'rs1283'
#galname = 'rf0503'
#galname = 'rf0376'
folderlist = {'rs0010': 'DS_2012-04-15/',
          'rs0346': 'EH_2014-09-20/',
          'rf0477': 'KE_2012-07-26/'}
folder = folderlist[galname]
folder = root+folder
filename = folder+'slincen_'+galname+'bspec.fits'
comp_file = folder+'cen_arc_broad.fits'
galaxy = fits.open(filename)
gal = galaxy[0].data
hgal = galaxy[0].header

#Wavelength solution of the galaxy
#NAXIS = Number of pixels on the x-axis (wavelength axis)
#CRPIX = Starting pixel number
#CDELT1 = delta(lambda) in Angstroms per pixel
#CRVAL1 = Wavelength value corresponding to starting pixel (CRPIX)
lam = (np.arange(hgal['NAXIS1']) + hgal['CRPIX1']-1) * hgal['CDELT1'] + hgal['CRVAL1'] 

standard_name = '/srv/one/resolve/analysis_code/emlineanalysis/broademline/R2013-11-04_LTT1788_cen_broad_sens.fits'
standard = fits.open(standard_name)

std = standard[0].data
hstd = standard[0].header

#Wavelength solution of the standard star
stdlam = (np.arange(hstd['NAXIS1']) + (hstd['CRPIX1']-1))* hstd['CDELT1'] + hstd['CRVAL1'] 

#Interpolate flux calibration to lam (galaxy spectrum wavelengths)
#The wavelength solution of the Standard Star and Galaxy may be differernt, 
#so resample the standard star flux to match the galaxy's wavelengths
f = interpolate.interp1d(stdlam,std, fill_value = 'extrapolate')
flux = f(lam)
fluxcorr = 10**(-0.4*(flux - np.median(flux))) #if sens in mag, gives ratios

#Multiply galaxy spectrum by fluxcorr for relative calibration (not absolute)
galcorr = gal*fluxcorr

#plt.figure()
#plt.plot(stdlam, std,'b')
#plt.plot(lam,flux,'r')
#plt.plot(lam,gal[int(len(gal)/2)],'g')
#plt.plot(lam,galcorr[int(len(gal)/2)],'k')


#Quick and dirty way to find the center
#center = np.where(gal[:, int(np.shape(gal)[1]/2)] == \
#                   max(gal[:, int(np.shape(gal)[1]/2)]))[0]

#If the number of rows in the galaxy is not a multiple of 6, zero pad rows to
#the closest multiple of 6.
#This is done so that the center of the new binned array is an integer number
shape = np.shape(galcorr)
if shape[0]%6:
    rows = 6 - (shape[0]%6)
    galcorr = np.append(galcorr,np.zeros([rows,shape[1]]), axis = 0)
shape = np.shape(galcorr)

#Finding center of the galaxy by finding mean of gaussian fit to the continuum
#center = findcenter(galcorr)
center = 180 #SDSS matched center
#Create an empty array to store the binned galaxy with 1/3rd the number of rows
#binned_gal = np.zeros([shape[0]//3, shape[1]])
#cen_bin = np.shape(binned_gal)[0]//2
#
##Binning from the center of the galaxy; every row in the binned galaxy gets the
##value of the total flux in 3 rows of the original flux array. This is done
##since there is no independent information in 3 adjacent rows of the original
##spectrum.
#for i in np.arange(0,np.shape(binned_gal)[0]//2):
#    binned_gal[cen_bin+i] = np.sum(galcorr[(center+(3*i) -1):(center+(3*i) +1)], axis = 0)
#
#    binned_gal[cen_bin-i] = np.sum(galcorr[(center-(3*i) -1):(center-(3*i) +1)], axis = 0)
#
#binned_gal = binned_gal[cen_bin-(binned_gal.shape[0]//5) : cen_bin+(binned_gal.shape[0]//5)]
binned_gal = np.zeros([1, shape[1]])
binned_gal[0,:] = np.sum(galcorr[(center-3):(center+3)],axis = 0)
#binned_gal = binned_gal[None,:]
#Add a dummy dimension and make the longslit spectrum a 3D data-cube
#Add keywords to the header corresponding to the 3rd dimension as wavelength
hdu = fits.getheader(filename)
binned_gal2 = binned_gal.T[:,None]
hdu['NAXIS'] = 3
hdu['NAXIS1'] = binned_gal2.shape[1]
hdu['NAXIS2'] = binned_gal2.shape[2]
hdu['NAXIS3'] = binned_gal2.shape[0]

hdu['CRVAL3'] = hdu['CRVAL1']
hdu['CRVAL1'] = 1

hdu['CRPIX3'] = 1

hdu['CDELT3'] = hdu['CDELT1']
hdu['CD3_3'] = hdu['CD1_1']
hdu['CDELT1'] = 1
hdu['CD1_1'] = 1

#outputfile = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/binned3d'+galname+'crop.fits'
outputfile = '/afs/cas.unc.edu/users/m/u/mugpol/Desktop/gistTutorial/inputData/sdss2arc'+galname+'.fits'

if not os.path.exists(outputfile):
    fits.writeto(outputfile, binned_gal2, hdu)


comp_file = folder+'lincen_'+galname+'bspec.fits'
compfile = comp_file
tmp1=comp_file.split('/')
tmp2=tmp1[7].split('.')
tmp3=tmp2[1]
grat=hgal['GRATING']
camang=hgal['CAM_ANG']
grtang=hgal['GRT_ANG']

lines = [5577.334, 6300.304, 6363.780, 6863.955]
#lines = [5578.5,5894.6,6301.7]#,7246.0]
#lines=[3650,4358,5460,6965]
waverange = hgal['CRVAL1'] + np.array([0.,hgal['CDELT1']*(hgal['NAXIS1']-1)])
add=[30, 30, 25, 10]
#add=[20, 20, 20, 20]

temp=fits.open(comp_file)
arcspec=temp[0].data
shap=np.shape(arcspec)

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
plt.plot(lam,arcspec[shap[0]//2+100,:],'b-')

plt.plot(lam[x1],arcspec[shap[0]//2+100,pixlines[0]-add[0]:pixlines[0]+add[0]],'r-')
plt.plot(lam[x2],arcspec[shap[0]//2+100,pixlines[1]-add[1]:pixlines[1]+add[1]],'r-')
plt.plot(lam[x3],arcspec[shap[0]//2+100,pixlines[2]-add[2]:pixlines[2]+add[2]],'r-')
plt.plot(lam[x4],arcspec[shap[0]//2+100,pixlines[3]-add[3]:pixlines[3]+add[3]],'r-')
#plt.plot(lam[x5],arcspec[shap[0]//2+100,pixlines[4]-add[4]:pixlines[4]+add[4]],'r-')


plt.title('red = arc lines used for FWHM calc')

fwhminpix=[]
for i in np.arange(len(lines)):
    nput=arcspec[shap[0]//2+100,pixlines[i]-add[i]:pixlines[i]+add[i]]
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






























#
#
#
#
#
#
## =====================================================================================
#
#def DER_SNR(flux):
#   
## =====================================================================================
#   from numpy import array, where, median, abs 
#
#   flux = array(flux)
#
#   # Values that are exactly zero (padded) are skipped
#   flux = array(flux[where(flux != 0.0)])
#   n    = len(flux)      
#
#   # For spectra shorter than this, no value can be returned
#   if (n>4):
#      signal = median(flux)
#
#      noise  = 0.6052697 * median(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
#      return noise #float(signal / noise)  
#
#   else:
#
#      return 0.0
#
## end DER_SNR -------------------------------------------------------------------------
#s     = np.shape(binned_gal2)
#print(s)
#spec  = np.reshape(binned_gal2,[s[0],s[1]*s[2]])
#espec = np.zeros(np.shape(spec))
#for i in range(spec.shape[1]):
#    espec[:,i] = DER_SNR(spec[:,i])
#wave = hgal['CRVAL1']+(np.arange(s[0]))*hgal['CD1_1']
#cvel      = 299792.458
#wave      = wave / (1+0.0219876)
#print(wave)
#idx       = np.where( np.logical_and( wave >= 3050, wave <= 7000 ) )[0]
#spec      = spec[idx,:]
#espec     = espec[idx,:]
#wave      = wave[idx]
#  
#idx_good = np.where( np.median(spec, axis=0) > 0.0 )[0]
#spec     = spec[:, idx_good]
#espec    = espec[:,idx_good]
#signal = np.nanmedian(spec,axis=0)
#noise = espec[0,:]
#snr    = signal / noise
#print(snr)