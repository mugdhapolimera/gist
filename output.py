# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 18:05:54 2019

@author: mugdhapolimera
"""

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.table import Table

folder = 'F:/mugdhapolimera/Shared VM/binned3drf0477crop_snr3test2/'
fit = Table.read(folder+'binned3drf0477crop_ppxf-bestfit.fits')
data = fits.open(folder+'binned3drf0477crop.fits')[0].data
hgal = fits.open(folder+'binned3drf0477crop.fits')[0].header
lam = (np.arange(hgal['NAXIS3']) + hgal['CRPIX3']-1) * hgal['CDELT3'] + hgal['CRVAL3'] 

plam = fits.open(folder+'binned3drf0477crop_ppxf-optimalTemplates.fits')[2].data['LOGLAM_TEMPLATE']
good = (lam > 4300) & (lam < 6900)
temp = fits.open(folder+'binned3drf0477crop_ppxf-optimalTemplates.fits')
plt.figure()
plt.plot(np.exp(plam), np.transpose(temp[1].data[0]))

#plt.plot(np.exp(plam), fit[0]['BESTFIT'])
plt.plot(lam[good],data[:,:,24][good]/max(data[:,:,24][good]))
#plt.plot(data[:,:,23][good])
#plt.plot(data[:,:,22][good])