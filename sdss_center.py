#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:24:21 2019

@author: mugpol

Find WCS solution for SOAR data and find the SDSS fibre center
"""

import astropy.io.fits as fits
import numpy as np
#import astropy.wcs as wcs
from scipy.io import readsav
import matplotlib.pyplot as plt
from scipy import spatial

folder = '/srv/two/resolve/working/spec/DS_2012-04-15/'
hdu = fits.open(folder+'slincen_rs0010bspec.fits')[0].header
#w = wcs.WCS(hdu)
#
#print(w.wcs.name)
#w.wcs.print_contents()
coords= readsav(folder+'coords_rs0010bspec.sav')
ra_cen = coords.coordra[np.where(coords.coordslice == 1)]
dec_cen = coords.coorddec[np.where(coords.coordslice == 1)]
x_cen = coords.coordx[np.where(coords.coordslice == 1)]
y_cen = coords.coordy[np.where(coords.coordslice == 1)]
wcs = list(zip(y_cen,ra_cen,dec_cen))
cont_cen = 181 #from findcenter function : continuum center
sdss_ra = 135.21266 #coordinates of the fiber
sdss_dec = 1.16109
print("Continuum Center", wcs[np.where(y_cen == cont_cen)[0][0]])

res_ra,res_dec = (135.212560, 1.161232)
tree = spatial.KDTree(list(zip(ra_cen,dec_cen)))
#closest = tree.query([(res_ra,res_dec)])
closest = tree.query([(sdss_ra,sdss_dec)])

#close_ra = min(ra_cen, key=lambda x:abs(x-sdss_ra))
#close_dec = min(dec_cen, key=lambda x:abs(x-sdss_dec))
#print(close_ra,close_dec)
#print(np.where(ra_cen == close_ra),np.where(dec_cen == close_dec))
print("SDSS Center", (sdss_ra,sdss_dec), closest)
closest = tree.query([(res_ra,res_dec)])
print("RESOLVE Center", (res_ra,res_dec), closest)

#fig,ax = plt.subplots()
#ax.plot(x_cen,y_cen,'o')
#
#for i,txt in enumerate(wcs):
#    ax.annotate(txt, (x_cen[i],y_cen[i]))