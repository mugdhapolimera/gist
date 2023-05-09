# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:24:09 2023

@author: mugdhapolimera
Compare decals and resolve r-band rad_50
"""
from scipy.optimize import curve_fit
import pandas as pd
from astropy.table import Table as t
from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.coordinates import * #SkyCoord
import matplotlib.pyplot as plt

def Gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

decals = t.read('decalsdr9_sdssdr17_xmatch.fits')
decals = decals.to_pandas()
decals = decals[(decals.shape_r > 0)]

fullres = pd.read_csv('RESOLVE_live2Oct2022.csv')
fullres.index = fullres.name

decalscoord = SkyCoord(ra=decals.ra*u.degree, dec=decals.dec*u.degree)
rescoord = SkyCoord(ra=fullres.radeg*u.degree, dec=fullres.dedeg*u.degree)
idx, d2d, d3d = rescoord.match_to_catalog_sky(decalscoord)

res_match = fullres[d2d < 5*u.arcsec]
res_match.index = np.arange(len(res_match))

idx = idx[d2d < 5*u.arcsec]
d3d = d3d[d2d < 5*u.arcsec]
d2d = d2d[d2d < 5*u.arcsec]
decalsres_match = decals.iloc[idx]
decalsres_match.index = np.arange(len(res_match))

decalsres = res_match.merge(decalsres_match,left_index = True, right_index = True)

radiuscol = 'shape_r'

r50good = (decalsres['pf_r50_arcseconds'] > 0) & (decalsres[radiuscol] > 0) & np.isfinite(decalsres[radiuscol])
percentdiff = (decalsres[radiuscol]- decalsres['pf_r50_arcseconds'])/decalsres['pf_r50_arcseconds']

ratio = decalsres[radiuscol][r50good]/decalsres['pf_r50_arcseconds'][r50good]

ratiodist = np.histogram(ratio, bins = 'fd')
bincen = (ratiodist[1][:-1] + ratiodist[1][1:])/2.0
parameters, covariance = curve_fit(Gauss, bincen , ratiodist[0])
norm = parameters[0]
mean = parameters[1]
sigma = parameters[2]
  
xaxis = np.arange(bincen[0], bincen[-1], step = abs(sigma)/4)
fit_y = Gauss(xaxis, norm, mean, sigma)

plt.figure()
plt.plot(xaxis, fit_y, color = 'k', lw = 3, zorder = 10)
plt.hist(ratio, bins = 'fd', color = 'b', histtype = 'step', \
         lw = 3, hatch = 'x')
