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

root = ''
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
        #pdb.set_trace()
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
lam, arcspec, conv, pix_to_fit = extract_spectra_parameters(skyfiles,arcfiles,'blue')

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

#plt.plot(lam,skyspec)
#plt.show()
