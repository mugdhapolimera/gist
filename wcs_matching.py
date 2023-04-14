# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 09:57:34 2019

@author: mugdhapolimera
"""
import numpy as np
import pandas as pd
from scipy.io import readsav
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os
#from astropy.wcs import WCS as wcs 
from astropy import wcs
from astropy.coordinates import SkyCoord as sky
from astropy import units as u
from matplotlib.collections import PatchCollection

resdata = readsav(r'C:\Users\mugdhapolimera\github\SDSS_spectra\resolvecatalog.dat')
galname = 'rs0010' #'rs0775'
if galname == 'rs0775':
    folder = '71146' #rs0775
    sdss= [183.25125,0.04375]
    resolve = [resdata['ra'][np.where(resdata['name'] == galname)][0], \
               resdata['dec'][np.where(resdata['name'] == galname)][0]]
    res_cont = [135.2128,1.1614498]

if galname == 'rs0010':
    folder = '372320' #rs0010
    sdss= [135.21266,1.16109]
    resolve = [resdata['ra'][np.where(resdata['name'] == galname)][0], \
               resdata['dec'][np.where(resdata['name'] == galname)][0]]
    res_cont = [135.2128,1.1614498]

v = resdata['vhel'][np.where(resdata['name'] == galname)] #km/s
z = v/3e5

os.chdir(r'F:\mugdhapolimera\Documents\UNC\Courses\Research\SAMI Data\/'+ folder)
image = folder+'_cube_red.fits'
cubehdu = fits.open(image)[0].header
lam0 = cubehdu['CRVAL3']-((cubehdu['CRPIX3']-1)*cubehdu['CDELT3'])
lam = (np.arange(cubehdu['NAXIS3']))*cubehdu['CDELT3'] + lam0
cubehdu['CDELT1'] = abs(cubehdu['CDELT1'])
cube = fits.open(image)[0].data#.flatten()


res = readsav(r'coords_'+galname+'bspec.sav')
cen = 1
ra_cen = res.coordra[np.where(res.coordslice == cen)]
dec_cen = res.coorddec[np.where(res.coordslice == cen)]
#sdss = sky(ra = 135.21266*u.degree,dec = 1.16109*u.degree, frame = 'icrs')
#sdsspix = w.wcs_world2pix([sdss],1)
lam_start = 6350*(1+z) 
lam_end = 6500*(1+z)
start = np.where(abs(lam-lam_start) == min(abs(lam-lam_start)))[0][0]
end = np.where(abs(lam-lam_end) == min(abs(lam-lam_end)))[0][0]
cube[np.where(np.isnan(cube))] = 0
w = wcs.WCS(cubehdu)
w = w.dropaxis(2)
xy = np.meshgrid(np.arange(50), np.arange(50))
coords = w.all_pix2world(xy[0], xy[1],0)
fig = plt.figure()
ax = plt.subplot(projection=w)
ax.imshow(np.sum(cube[start:end,:,:], axis = 0), origin='lower', cmap=plt.cm.viridis)
#ax.scatter(sdsspix[0][0], sdsspix[0][1],color = 'm')
resC = [Circle((ra,dec), radius = 0.5/3600.0) for ra,dec in zip(ra_cen,dec_cen)]
#respatches = PatchCollection(resC, edgecolor = 'r', facecolor = 'none',
#                             transform = ax.get_transform('world'), 
#                             label = 'RESOLVE WCS Coordinates')
c = Circle((sdss[0],sdss[1]), radius = 1.5/3600.0, edgecolor = 'magenta', lw = 5,
           facecolor = 'none', transform = ax.get_transform('world'), 
           label = 'SDSS 2" fibre')
ax.add_patch(c)
#ax.add_collection(respatches)
#ax.scatter(sdss[0],sdss[1], c = 'm', s = 100,
#           transform = ax.get_transform('world'), label = 'SDSS Center')
#ax.scatter(resolve[0],resolve[1], c = 'r', s = 100,
#           transform = ax.get_transform('world'), 
#           label = 'RESOLVE Photometric Center')
#ax.scatter(res_cont[0],res_cont[1], c = 'w', marker = '^', s = 100, 
#           edgecolors = 'k', linewidths = 3, transform = ax.get_transform('world'), 
#           label = 'RESOLVE Continuum Center')
plt.xlabel('RA', fontsize = 22)
plt.ylabel('Dec', fontsize = 22)
ra_low = ra_cen - 0.5/3600.0
ra_high = ra_cen + 0.5/3600.0
dec_low = dec_cen - 0.5/3600.0
dec_high = dec_cen + 0.5/3600.0

slitlow = np.poly1d(np.polyfit(ra_low,dec_low,30))
slithigh = np.poly1d(np.polyfit(ra_high,dec_high,30))
ra = np.arange(ra_low.min(), ra_low.max(),0.0001)
dec = slitlow(ra)
#ax.plot(ra,dec,c = 'green', lw = 5,
#           transform = ax.get_transform('world'))
#ra = np.arange(ra_high.min(), ra_high.max(),0.0001)
#dec = slithigh(ra)
#ax.plot(ra,dec,c = 'green', lw = 5,
#           transform = ax.get_transform('world'), label = 'SOAR 1" Slit')
pix = np.where((coords[1] <= slitlow(coords[0])) & (coords[1] >= slithigh(coords[0])))
#ax.scatter(coords[0][pix],coords[1][pix], c = 'orange', s = 2,
#           transform = ax.get_transform('world'))
#plt.legend()

#plt.figure()
#plt.plot(ra_low, dec_low)
#plt.plot(ra,dec,'ro')
