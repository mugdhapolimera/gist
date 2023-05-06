# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 08:20:54 2021

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import sys
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from scipy.optimize import curve_fit
from bpt_s06_agn import s06_bptagn
from optimized_bpt_plotagn import bpt_plots
import matplotlib.gridspec as gridspec
os.chdir("../SDSS_spectra/mid_ir")
from midiragnplot import midir_plotallagn

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 


#Read in RESOLVE/ECO extinction corrected and S/N filtered data
save = 1
resolve = 0
eco = 0
full = 0
sdss = 1
#sdsscat = 'port'
sdsscat = 'jhu'
#sdsscat = 'nsa'
#sdsscat = 'master'
if sys.platform == 'linux2':
        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')

else:
    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')


if full: 
    if sdsscat == 'jhu':
        s06inputfile = 'ECO+RESOLVE_bpt1snr5_dext_jhu.csv'#'ECO_full_snr5.csv'
        bptinputfile = 'ECO+RESOLVE_bpt1snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
        selinputfile = 'ECO+RESOLVE_snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
    print 'ECO+RESOLVE RESULTS'
    s06outputfile = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
    bptoutputfile = 'eco+resolve_emlineclass_bpt1snr5_'+sdsscat+'.csv'
    midirfile = 'mid_ir/ECO+RESOLVE_WISE_good_final.csv'
    survey = 'ECO+RESOLVE'
elif eco: 
    if sdsscat == 'jhu':
        s06inputfile = 'ECO_full_bpt1snr5_dext_jhu.csv'#'ECO_full_snr5.csv'
        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
    if sdsscat == 'port':
        s06inputfile = 'ECO_full_bpt1snr5_dext_port.csv'#'ECO_full_snr5_port.csv'
        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_port.csv' #'ECO_full_bary_port.csv'#'ECO_full_snr5_port.csv'
    if sdsscat == 'nsa':       
        s06inputfile = 'ECO_full_bpt1snr5_dext_nsa.csv'
        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_nsa.csv'
    print 'ECO RESULTS'
    s06outputfile = 'eco_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
    bptoutputfile = 'eco_emlineclass_dext_snr5_'+sdsscat+'.csv'
    inobssamplefile = 'ECO_inobssample.csv'
    midirfile = 'mid_ir/ECO_WISE_good.csv'
    survey = 'ECO'

elif resolve:
    if sdsscat == 'port':       
        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_port.csv'
        bptinputfile = 'RESOLVE_full_snr5_dext_port.csv'
    if sdsscat == 'jhu':       
        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_jhu.csv'
        bptinputfile = 'RESOLVE_full_snr5_dext_jhu.csv'
    if sdsscat == 'nsa':       
        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_nsa.csv'
        bptinputfile = 'RESOLVE_full_snr5_dext_nsa.csv'
    if sdsscat == 'master':       
        bptinputfile = 'RESOLVE_snr5_master_new.csv'
        s06inputfile = 'RESOLVE_snr5_master_new.csv'
    print 'RESOLVE RESULTS'
    s06outputfile = 'resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
    bptoutputfile = 'resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
    inobssamplefile = 'RESOLVE_inobssample.csv'
    survey = 'RESOLVE'
#    midirfile = 'mid_ir/RESOLVE_WISE_good_randerr.csv'
    midirfile = 'mid_ir/ECO+RESOLVE_WISE_good_final.csv'

if sdss: 
    if sdsscat == 'jhu':
        s06inputfile = 'JHU-SDSSdr8_bpt1snr5_dext_jhu_new.csv'
        bptinputfile = 'JHU-SDSSdr8_bpt1snr5_dext_jhu_new.csv'
        selinputfile = 'JHU-SDSSdr8_snr5_dext_jhu.csv'
    midirfile = ''
    survey = 'JHU-SDSSdr8'
    s06outputfile = survey+'_s06emlineclass_dext_bpt1snr5_'+sdsscat+'.csv'
    bptoutputfile = survey+'_emlineclass_bpt1snr5_'+sdsscat+'.csv'
    print survey+' RESULTS'

#plt.figure()
#gs0 = gridspec.GridSpec(2, 1)
#gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0])
#gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1])
#
#ax1 = plt.subplot(gs00[0, 0])
#ax2 = plt.subplot(gs00[0, 1])

#ax3 = plt.subplot(gs01[0, 0])
#ax4 = plt.subplot(gs01[0, 1])
#ax5 = plt.subplot(gs01[0, 2])
#fig, ax1 = plt.subplots(1,1, sharey = True)
#
#ax1 = s06_bptagn(s06inputfile, s06outputfile, bptoutputfile, midirfile,
#              eco, resolve, full, sdsscat, save, ax1)

#midirfile = 'mid_ir/GAMA_WISECat.fits'
#fig, ax2 = plt.subplots(1,1, sharey = True)
#ax2 = midir_plotallagn(ax2,midirfile, s06outputfile, bptoutputfile, survey, save)
#os.chdir("C:\Users\mugdhapolimera\github\SDSS_Spectra\/")
#

simple = 1
fig, (ax3, ax4, ax5) = plt.subplots(1,3, sharey = True)

ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile,  selinputfile, s06outputfile,  midirfile,
                          eco, resolve, full, sdsscat, save, ax3, ax4, ax5, simple)

#plt.show()

#=======
## -*- coding: utf-8 -*-
#"""
#Created on Tue Mar 22 13:00:52 2022
#
#@author: mugdhapolimera
#"""
#
## -*- coding: utf-8 -*-
#"""
#Created on Fri Jun 11 08:20:54 2021
#
#@author: mugdhapolimera
#"""
#
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import os 
#import sys
#from scipy.stats import kde
#import matplotlib
#import matplotlib.cm as cm
#import matplotlib.colors as colors
#matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#from scipy.optimize import curve_fit
#from bpt_s06_agn import s06_bptagn
#from optimized_bpt_plotagn import bpt_plots
#import matplotlib.gridspec as gridspec
#os.chdir("../SDSS_spectra/mid_ir")
#from midiragnplot import midir_plotallagn
#
#def density_estimation(m1, m2):
#    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
#    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
#    values = np.vstack([m1, m2])                                                                        
#    kernel = scipy.stats.gaussian_kde(values)                                                                 
#    Z = np.reshape(kernel(positions).T, X.shape)
#    return X, Y, Z 
#
#
##Read in RESOLVE/ECO extinction corrected and S/N filtered data
#save = 1
#resolve = 0
#eco = 0
#full = 1
##sdsscat = 'port'
#sdsscat = 'jhu'
##sdsscat = 'nsa'
##sdsscat = 'master'
#if sys.platform == 'linux2':
#        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')
#
#else:
#    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
#
#
#if full: 
#    if sdsscat == 'jhu':
#        s06inputfile = 'ECO+RESOLVE_bpt1snr5_dext_jhu.csv'#'ECO_full_snr5.csv'
#        bptinputfile = 'ECO+RESOLVE_bpt1snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
#        selinputfile = 'ECO+RESOLVE_snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
#    print 'ECO+RESOLVE RESULTS'
#    s06outputfile = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#    bptoutputfile = 'eco+resolve_emlineclass_dext_updatedsample_'+sdsscat+'.csv'
#    midirfile = 'mid_ir/RESOLVE_WISE_good_randerr.csv'
#    survey = 'ECO+RESOLVE'
#elif eco: 
#    if sdsscat == 'jhu':
#        s06inputfile = 'ECO_full_bpt1snr5_dext_jhu.csv'#'ECO_full_snr5.csv'
#        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_jhu.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
#    if sdsscat == 'port':
#        s06inputfile = 'ECO_full_bpt1snr5_dext_port.csv'#'ECO_full_snr5_port.csv'
#        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_port.csv' #'ECO_full_bary_port.csv'#'ECO_full_snr5_port.csv'
#    if sdsscat == 'nsa':       
#        s06inputfile = 'ECO_full_bpt1snr5_dext_nsa.csv'
#        bptinputfile = 'ECO/SEL/ECO_full_snr5_dext_nsa.csv'
#    print 'ECO RESULTS'
#    s06outputfile = 'eco_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#    bptoutputfile = 'eco_emlineclass_dext_snr5_'+sdsscat+'.csv'
#    inobssamplefile = 'ECO_inobssample.csv'
#    midirfile = 'mid_ir/ECO_WISE_good.csv'
#    survey = 'ECO'
#
#else:
#    if sdsscat == 'port':       
#        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_port.csv'
#        bptinputfile = 'RESOLVE_full_snr5_dext_port.csv'
#    if sdsscat == 'jhu':       
#        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_jhu.csv'
#        bptinputfile = 'RESOLVE_full_snr5_dext_jhu.csv'
#    if sdsscat == 'nsa':       
#        s06inputfile = 'RESOLVE_full_bpt1snr5_dext_nsa.csv'
#        bptinputfile = 'RESOLVE_full_snr5_dext_nsa.csv'
#    if sdsscat == 'master':       
#        bptinputfile = 'RESOLVE_snr5_master_new.csv'
#        s06inputfile = 'RESOLVE_snr5_master_new.csv'
#    print 'RESOLVE RESULTS'
#    s06outputfile = 'resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#    bptoutputfile = 'resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
#    inobssamplefile = 'RESOLVE_inobssample.csv'
#    survey = 'RESOLVE'
##    midirfile = 'mid_ir/RESOLVE_WISE_good_randerr.csv'
#    midirfile = 'mid_ir/RESOLVE_WISE_good_syserr.csv'
#
##plt.figure()
##gs0 = gridspec.GridSpec(2, 1)
##gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0])
##gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1])
##
##ax1 = plt.subplot(gs00[0, 0])
##ax2 = plt.subplot(gs00[0, 1])
#
##ax3 = plt.subplot(gs01[0, 0])
##ax4 = plt.subplot(gs01[0, 1])
##ax5 = plt.subplot(gs01[0, 2])
#fig, ax1 = plt.subplots(1,1, sharey = True)
#
#ax1 = s06_bptagn(s06inputfile, s06outputfile, bptoutputfile, midirfile,
#              eco, resolve, full, sdsscat, save, ax1)
#
##midirfile = 'mid_ir/GAMA_WISECat.fits'
##fig, ax2 = plt.subplots(1,1, sharey = True)
##ax2 = midir_plotallagn(ax2,midirfile, s06outputfile, bptoutputfile, survey, save)
##os.chdir("C:\Users\mugdhapolimera\github\SDSS_Spectra\/")
##
#simple = 1
#fig, (ax3, ax4, ax5) = plt.subplots(1,3, sharey = True)
#
#ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile,  selinputfile, s06outputfile,  midirfile,
#                          eco, resolve, full, sdsscat, save, ax3, ax4, ax5, simple, allcat = 0)
#
##plt.show()
#
#>>>>>>> 4d300bdd936fd05600c51bbe79ca4a5a518e4724
