# -*- coding: utf-8 -*-
"""
Creating a new search catalog for nuggets and AGN in the local universe
Decals photometry - color, radius
SDSS spectroscopy - MPA-JHU
GSLWC - SFRs

Created on Thu Mar 30 13:52:49 2023

@author: mugdhapolimera
"""

import pandas as pd
from astropy.table import Table as t
from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.coordinates import * #SkyCoord
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP5

#for decals API
#from getpass import getpass

crossmatchedfile = 'decalsdr7-JHU-GSWLC-crossmatch.csv'
crossmatch = 0

if crossmatch:
    #DECALS MATCHED WITH SDSS DR14
#    decals = pd.read_csv('decals_dr7_sdssdr14.csv')
#    decals = decals[decals.decals_id > 0]
    
    decals = t.read('decalsdr9_sdssdr17_xmatch.fits')
    decals = decals.to_pandas()
    decals = decals[decals.shape_r > 0]

    
    #MPA-JHU EMISSION LINE MEASUREMENTS
    jhu = t.read('galSpecLine-dr8.fits')
    jhu = jhu.to_pandas()
    jhuinfo = t.read('galSpecInfo-dr8.fits')
    jhuextra = t.read('galSpecExtra-dr8.fits')
    sdss = fits.open('sdssdr8_phot_redshift_mugdhapolimera.fit')[1].data
    sdss = pd.DataFrame(sdss)
    for x in sdss.keys():
        sdss[x] = np.array(sdss[x]).byteswap().newbyteorder()
    sdss['SPECOBJID'] = np.array(sdss['specObjID'], dtype = np.float64)
    
    jhu['logmstar'] = np.array(jhuextra['LGM_TOT_P50'], '<f8')
    jhu['jhura'] = np.array(jhuinfo['RA'], '<f8')
    jhu['jhudec'] = np.array(jhuinfo['DEC'], '<f8')
    jhu['reliable'] = np.array(jhuinfo['RELIABLE'], '<f8')
    jhu['targettype'] = np.array(jhuinfo['TARGETTYPE'], 'str')
    jhu['jhuredshift'] = np.array(jhuinfo['Z'], '<f8')

    jhu = jhu[jhu.SPECOBJID != '                   ']
    jhu = jhu.astype({'SPECOBJID': np.float64})        
    jhu = jhu[jhu.jhura != -9999]
    jhu = jhu[jhu.reliable == 1] 
    jhu = jhu[jhu.targettype.str.strip() == 'GALAXY']
    jhu = jhu[jhu.logmstar > 0]

    jhu = jhu.merge(sdss, how = 'inner', on = 'SPECOBJID')
    jhusdsstab= t.from_pandas(jhu)
    jhusdsstab.write('JHU-SDSSdr8-xmatch.fits')


    #GSWLC CATALOG WITH SFRs
    salimcols = [u'objid', u'glxid', u'plate',
           u'mjd', u'fiberid', u'ra', u'dec', u'z', u'chisq_r', u'salim_logmstar',
           u'err_salim_logmstar', u'logsfr_sed', u'err_logsfr_sed', u'a_fuv',
           u'err_a_fuv', u'a_b', u'err_a_b', u'a_v', u'err_a_v', u'flag_sed',
           u'uvsurvey', u'logsfr_ir_wise', u'flag_wise', u'logsfr_ir_unwise',
           u'flag_unwise', u'flag_mgs']
    salim = pd.read_csv('GSWLC-X1.dat', sep = '\s+', names = salimcols)
    salim = salim[salim.logsfr_sed > -99]
    
    
    #DECALS and MPA-JHU CROSSMATCH
    decalscoord = SkyCoord(ra=decals.ra*u.degree, dec=decals.dec*u.degree)
    jhucoord = SkyCoord(ra=jhu.jhura*u.degree, dec=jhu.jhudec*u.degree)
    idx, d2d, d3d = jhucoord.match_to_catalog_sky(decalscoord)
    
    jhu_match = jhu[d2d < 1.5*u.arcsec]
    jhu_match.index = np.arange(len(jhu_match))
    
    idx = idx[d2d < 1.5*u.arcsec]
    d3d = d3d[d2d < 1.5*u.arcsec]
    d2d = d2d[d2d < 1.5*u.arcsec]
    decals_match = decals.iloc[idx]
    decals_match.index = np.arange(len(decals_match))
    
    decals_jhu = decals_match.merge(jhu_match,left_index = True, right_index = True)
#    decals_jhu['ra'] = decals_jhu['ra_y']
#    decals_jhu['dec'] = decals_jhu['dec_y']
#    decals_jhu['ra_decals'] = decals_jhu['ra_x']
#    decals_jhu['dec_decals'] = decals_jhu['dec_x']
    
    mapping = {'ra_y': 'ra', 'ra_x': 'ra_decals','dec_y': 'dec', 'dec_x': 'dec_decals'}
    decals_jhu = decals_jhu.rename(columns=mapping)
    
    decals_jhutab = t.from_pandas(decals_jhu)
    decals_jhutab.write('decalsdr9-SDSSdr8-xmatch.fits')
    
    #DECALS+MPA-JHU and GSWLC CROSSMATCH
    decalsjhucoord = SkyCoord(ra=decals_jhu.jhura*u.degree, dec=decals_jhu.jhudec*u.degree)
    salimcoord = SkyCoord(ra=salim.ra*u.degree, dec=salim.dec*u.degree)
    idx, d2d, d3d = salimcoord.match_to_catalog_sky(decalsjhucoord)
    
#    idx, d2d, d3d = decalsjhucoord.match_to_catalog_sky(salimcoord)
#    decalsjhu_match = decals_jhu[d2d < 1.5*u.arcsec]
#    decalsjhu_match.index = np.arange(len(decalsjhu_match))

    salim_match = salim[d2d < 1.5*u.arcsec]
    salim_match.index = np.arange(len(salim_match))
    
    idx = idx[d2d < 1.5*u.arcsec]
    d3d = d3d[d2d < 1.5*u.arcsec]
    d2d = d2d[d2d < 1.5*u.arcsec]
    decalsjhu_match = decals_jhu.iloc[idx]
    decalsjhu_match.index = np.arange(len(decalsjhu_match))
#    salim_match = salim.iloc[idx]
#    salim_match.index = np.arange(len(salim_match))
    
    decalsjhusalim = decalsjhu_match.merge(salim_match,left_index = True, right_index = True)
    mapping = {'ra_y': 'ra_salim', 'ra_x': 'ra','dec_y': 'dec_salim', 'dec_x': 'dec'}
    decalsjhusalim = decalsjhusalim.rename(columns=mapping)
    
    decalsjhusalimtab = t.from_pandas(decalsjhusalim)
    decalsjhusalimtab.write('decalsdr9-SDSSdr8-GSWLC-xmatch.fits')

decalsjhusalim = t.read('decalsdr9-SDSSdr8-GSWLC-xmatch.fits')
decalsjhusalim= decalsjhusalim.to_pandas()
decalsjhusalim = decalsjhusalim[(decalsjhusalim.jhuredshift > 0) & (decalsjhusalim.r > 0)]
cosmo = FlatLambdaCDM(H0 = 70*u.km/u.s/u.Mpc, Om0=0.3)
decalsjhusalim['distance'] = cosmo.luminosity_distance(decalsjhusalim['jhuredshift']) #in Mpc
distmod = np.array(Distance(decalsjhusalim['distance'], unit= u.Mpc).distmod)
decalsjhusalim['absrmag'] = decalsjhusalim['r'] - distmod
#decalsjhusalim = pd.read_csv('decalsdr7-JHU-GSWLC-crossmatch.csv')

#Filter good:
#good masses and radii
#search = decalsjhusalim[decalsjhusalim.logmstar > 0]
jhu = t.read('JHU-SDSSdr8-xmatch.fits')
jhu = jhu.to_pandas()
zcol = 'jhuredshift'
jhu = jhu[(jhu[zcol]> 0) & (jhu.r > 0)]
#
#jhu = jhu[(jhu[zcol]> 0.005) & (jhu[zcol] < 0.066) & (jhu.r > 0)]

###############################################################################
# Defining Volume-limited sample
###############################################################################

rlim = 17.77
M_rmax = -18
z_uplim = 10**((rlim - M_rmax + 5)/5.0)/1e6 * 70/3e5

M_rmax2 = -19.5
z_uplim2 = 10**((rlim - M_rmax2 + 5)/5.0)/1e6 * 70/3e5

cosmo = FlatLambdaCDM(H0 = 70*u.km/u.s/u.Mpc, Om0=0.3)
jhu['distance'] = cosmo.luminosity_distance(jhu[zcol]) #in Mpc
distmod = np.array(Distance(jhu['distance'], unit= u.Mpc).distmod)
jhu['absrmag'] = jhu['r'] - distmod


#plt.figure()
#plt.plot(jhu[zcol],jhu['absrmag'],',')
#plt.ylim(-17,-24)
#plt.axhline(M_rmax,color = 'k')
#plt.axvline(0.005, color = 'k')
#plt.axvline(z_uplim, color = 'k')
#plt.xlabel('z')
#plt.ylabel('M_r')
#
#plt.figure()
#plt.plot(jhu[zcol],jhu['absrmag'],',')
#plt.axhline(M_rmax2,color = 'r', ls = '--')
#plt.axvline(0.005, color = 'r')
#plt.axvline(z_uplim2, color = 'r')
#plt.xlabel('z')
#plt.ylabel('M_r')
#plt.ylim(-17,-24)
#plt.axvspan(0.005,0.066, alpha = 0.15, color = 'k')
#plt.xlim(0,0.34)

jhu_zlim = jhu[(jhu[zcol] > 0.01) & (jhu[zcol] < 0.066)]
jhu_zlim = jhu_zlim.drop('SPECOBJID', axis = 1)
jhu_zlim.columns = jhu_zlim.columns.str.lower()
#jhu_zlim[['ra', 'dec', 'jhuredshift', 'absrmag']].to_csv('JHU-SDSSdr8-zlim-coords.csv')
#jhu_zlim.to_csv('JHU-SDSSdr8-zlim.csv')

decalsjhusalim.index = decalsjhusalim.SPECOBJID
search = decalsjhusalim.loc[jhu_zlim.specobjid]
search = search[(search.logmstar > 0) & (search.shape_r > 0)]

#optical AGN emission lines for BPT plot
searchopt = search[(search['H_ALPHA_FLUX'] > 0) & (search['H_ALPHA_FLUX'] > 5*search['H_ALPHA_FLUX_ERR']) & \
                   (search['H_BETA_FLUX'] > 0) & (search['H_BETA_FLUX'] > 5*search['H_BETA_FLUX_ERR']) & \
                   (search['OIII_5007_FLUX'] > 0) & (search['OIII_5007_FLUX'] > 5*search['OIII_5007_FLUX_ERR'])]

#mid-IR AGN using color
search['flux_w1_err'] = np.sqrt(1/search['flux_ivar_w1'])
search['flux_w2_err'] = np.sqrt(1/search['flux_ivar_w2'])
search['flux_w3_err'] = np.sqrt(1/search['flux_ivar_w3'])
searchir = search[(search['flux_w1'] > 0) & (search['flux_w1'] > 10*search['flux_w1_err']) & \
                  (search['flux_w2'] > 0) & (search['flux_w2'] > 10*search['flux_w2_err'])]

searchtab = t.from_pandas(search)
#searchtab.write('decalsdr9-SDSSdr8-nuggetsearchsample.fits', overwrite = True)
###############################################################################
print 'Nugget search sample: ', len(search)
print 'Optical AGN subsample: ', len(searchopt)
print 'IR AGN subsample: ', len(searchir)


#RESOLVE stats

#nuggets = pd.read_csv('resolveprofitcolorandsfrnuggets.csv')
#nuggets.index = nuggets.galname
#resolve = pd.read_csv('RESOLVE_barysample.csv')
#resolve.index = resolve.name
##resolve = resolve[resolve.logmstar < 9.5]
#nuggets = nuggets.loc[resolve.name]
#
#colnugget = (nuggets.bluenugget == 1) | (nuggets.greennugget == 1)| (nuggets.rednugget== 1)
#sfnugget = (nuggets.LSFnugget == 1)| (nuggets.MSFnugget == 1)| (nuggets.HSFnugget == 1)
#
#colnugname = set(nuggets.index[colnugget])
#sfnugname = set(nuggets.index[sfnugget])
#
#print '\nRESOLVE\n'
#print 'Color-based nuggets: ', np.sum(colnugget)
#print 'SF-based nuggets: ', np.sum(sfnugget)
#print 'common nuggets: ', np.sum(sfnugget & colnugget)
#
#print 'Optical AGN: ', len(optagn & set(resolve.name))
#print 'IR AGN: ', len(iragn & set(resolve.name))
#
#print 'Optical AGN in color-based nuggets: ', len(optagn & colnugname)
#print 'Optical AGN in SF-based nuggets: ', len(optagn & colnugname)
#
#print 'IR AGN in color-based nuggets: ', len(iragn & colnugname)
#print 'IR AGN in SF-based nuggets: ', len(iragn & colnugname)
#
#
###############################################################################
#DECALS and RESOLVE crossmatch
#decals = decalsjhusalim 
#decals = t.read('decalsdr9_sdssdr17_xmatch.fits')
#decals = decals.to_pandas()
#decals = decals[(decals.shape_r > 0)]
#
##nuggets = pd.read_csv('resolveprofitcolorandsfrnuggets.csv')
##nuggets.index = nuggets.galname
##fullres = pd.read_csv('RESOLVE_live02Feb2022.csv')
#fullres = pd.read_csv('ECO_live2Oct2022.csv')
#fullres.index = fullres.name
#
##nuggets = nuggets.loc[fullres.index]
##fullres['pf_r50_arcseconds'] = np.zeros(len(fullres))
##fullres['pf_r50_arcseconds'] = nuggets['pf_r50_arcseconds']
#
#decalscoord = SkyCoord(ra=decals.ra*u.degree, dec=decals.dec*u.degree)
#rescoord = SkyCoord(ra=fullres.radeg*u.degree, dec=fullres.dedeg*u.degree)
#idx, d2d, d3d = rescoord.match_to_catalog_sky(decalscoord)
#
#res_match = fullres[d2d < 5*u.arcsec]
#res_match.index = np.arange(len(res_match))
#
#idx = idx[d2d < 5*u.arcsec]
#d3d = d3d[d2d < 5*u.arcsec]
#d2d = d2d[d2d < 5*u.arcsec]
#decalsres_match = decals.iloc[idx]
#decalsres_match.index = np.arange(len(res_match))
#
#decalsres = res_match.merge(decalsres_match,left_index = True, right_index = True)
#
#devgood = (decalsres['shapedev_r'] > 0)
#plt.figure()
#plt.plot(decalsres['shapedev_r'][devgood], decalsres['pf_r50_arcseconds'][devgood], 'o')
#plt.xlabel('decals half-light radius from deVaucouleurs model (>0)')
#plt.ylabel('resolve profit r50')
#plt.plot(np.arange(160), np.arange(160), 'k')
#
#expgood = (decalsres['shapeexp_r'] > 0)
#plt.figure()
#plt.plot(decalsres['shapeexp_r'][expgood], decalsres['pf_r50_arcseconds'][expgood], 'o')
#plt.xlabel('decals half-light radius from exponential model (>0)')
#plt.ylabel('resolve profit r50')
#plt.plot(np.arange(160), np.arange(160), 'k')
#
##percentdiff = (decalsres['shapeexp_r'][expgood] - decalsres['pf_r50_arcseconds'][expgood])/decalsres['pf_r50_arcseconds'][expgood]
#percentdiff = decalsres['pf_r50_arcseconds'][expgood]/decalsres['shapeexp_r'][expgood]
#plt.figure()
#plt.plot(decalsres['pf_r50_arcseconds'][expgood], percentdiff,'o')
#plt.ylabel('(decals exponential r50 - resolve r50)/resolve r50')
#plt.xlabel('resolve profit r50') 
#plt.axhline(0,color = 'k')
#
#decalsres['decals_r50'] = np.zeros(len(decalsres))
#sel = np.sqrt(1/decalsres.shapedev_r_ivar) < np.sqrt(1/decalsres.shapeexp_r_ivar)
#decalsres['decals_r50'][sel] = decalsres['shapedev_r'][sel]
#decalsres['decals_r50'][~sel] = decalsres['shapeexp_r'][~sel]
#sel = decalsres.type == 'EXP'
#decalsres['decals_r50'][sel] = decalsres['shapeexp_r'][sel]
#sel = decalsres.type == 'DEV'
#decalsres['decals_r50'][sel] = decalsres['shapedev_r'][sel]

#radiuscol = 'shape_r'
##radiuscol = 'decals_r50'
#
#r50good = (decalsres['pf_r50_arcseconds'] > 0) & (decalsres[radiuscol] > 0) & np.isfinite(decalsres[radiuscol])
#percentdiff = (decalsres[radiuscol]- decalsres['pf_r50_arcseconds'])/decalsres['pf_r50_arcseconds']
##percentdiff = decalsres['pf_r50_arcseconds'][expgood]/decalsres['shapeexp_r'][expgood]
#plt.figure()
#plt.plot(decalsres['pf_r50_arcseconds'][r50good], percentdiff[r50good],'o')
#plt.ylabel('(decals r50 - res r50)/res r50')
#plt.xlabel('resolve profit r50') 
#plt.axhline(0,color = 'k')
#
#ratio = decalsres[radiuscol][r50good]/decalsres['pf_r50_arcseconds'][r50good]
#plt.figure()
#plt.plot(decalsres['pf_r50_arcseconds'][r50good], ratio, 'o')
#plt.ylabel('decals r50/ res r50')
#plt.xlabel('resolve profit r50') 
#plt.axhline(1,color = 'k')
#
#from scipy.optimize import curve_fit
#
#def Gauss(x, a, x0, sigma):
#    return a*np.exp(-(x-x0)**2/(2*sigma**2))
#
#ratiodist = np.histogram(ratio, bins = 'fd')
#bincen = (ratiodist[1][:-1] + ratiodist[1][1:])/2.0
#parameters, covariance = curve_fit(Gauss, bincen , ratiodist[0])
#norm = parameters[0]
#mean = parameters[1]
#sigma = parameters[2]
#  
#xaxis = np.arange(bincen[0], bincen[-1], step = abs(sigma)/4)
#fit_y = Gauss(xaxis, norm, mean, sigma)
#plt.figure()
#plt.plot(xaxis, fit_y, color = 'k', lw = 3, zorder = 10)
#plt.hist(ratio, bins = 'fd', color = 'b', histtype = 'step', \
#         lw = 3, hatch = 'x')
#
#plt.figure()
#plt.plot(decalsres[radiuscol][r50good], decalsres['pf_r50_arcseconds'][r50good], 'o')
#plt.ylabel('decals preferred r50')
#plt.xlabel('resolve profit r50') 
#plt.plot(np.arange(0,50), np.arange(0,50),color = 'k')
#
#
#
#plt.figure()
#plt.plot(decals['ra'], decals['dec_x'], 'o')
#plt.plot(decalsres['ra'], decalsres['dec_x'], 'o')
#plt.plot(fullres['radeg'], fullres['dedeg'], 'o', mfc = 'none', mec = 'k' )