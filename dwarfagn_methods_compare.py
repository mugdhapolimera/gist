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
from bpt_s06 import s06_bpt
from bpt_plots_new import bpt_plots
import matplotlib.gridspec as gridspec
os.chdir("../SDSS_spectra/mid_ir")
from midir_agn import midiragnplot
pd.set_option('display.max_columns', None)
def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 

#Read in RESOLVE/ECO extinction corrected and S/N filtered data
he2_flag = 0
save = 1
resolve = 0
eco = 0
full = 0
sdss = 1
#sdsscat = 'port'

sdsscat = 'jhu'
# = 'nsa'
#sdsscat = 'master'
if sys.platform == 'linux2':
        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')

else:
    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

namecol = 'name'
if full: 
    s06inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_bpt1snr5_dext_'+sdsscat+'.csv'#'ECO_full_snr5.csv'
    bptinputfile = 'ECO+RESOLVE_snr5_dext_'+sdsscat+'.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
    print 'ECO+RESOLVE RESULTS'
    s06outputfile = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
    bptoutputfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
    inobssamplefile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_barysample.csv'
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
    survey = 'ECO'

else:
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
if sdss: 
    if sdsscat == 'jhu':
        s06inputfile = 'JHU-SDSSdr8_bpt1snr5_dext_jhu_new.csv'
        bptinputfile = 'JHU-SDSSdr8_bpt1snr5_dext_jhu_new.csv'
        selinputfile = 'JHU-SDSSdr8_snr5_dext_jhu.csv'
    midirfile = ''
    survey = 'JHU-SDSSdr8'
    s06outputfile = survey+'_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
    bptoutputfile = survey+'_emlineclass_bpt1snr5_'+sdsscat+'.csv'
    print survey+' RESULTS'
    namecol = 'jhu_index'

plt.figure()
#gs0 = gridspec.GridSpec(2, 1)
#gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0])
#gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1])
#
#ax1 = plt.subplot(gs00[0, 0])
#ax2 = plt.subplot(gs00[0, 1])
#
#ax3 = plt.subplot(gs01[0, 0])
#ax4 = plt.subplot(gs01[0, 1])
#ax5 = plt.subplot(gs01[0, 2])
ax1 = plt.subplot(111)

ax1 = s06_bpt(s06inputfile, s06outputfile, eco, resolve, full, sdsscat, save, ax1, namecol)
#ax2 = plt.subplot(111)
#ax2 = midiragnplot(ax2,inobssamplefile, survey, save)

os.chdir("C:\Users\mugdhapolimera\github\SDSS_Spectra\/")

fig, (ax3, ax4, ax5) = plt.subplots(1,3, sharey = True)

ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile, eco, resolve, full, sdsscat, 
                          save, ax3, ax4, ax5)
#
#plt.show()

#from astropy.stats import binom_conf_interval
#pc = 100.0*k/n
#err_down, err_up = 100.0*binom_conf_interval(k, n) - pc

#inobssamplefile = bptinputfile
#plt.figure()    
#fig, (ax3, ax4, ax5) = plt.subplots(1,3, sharey = True)
#ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile, eco, resolve, full, sdsscat, 
#                          save, ax3, ax4, ax5)
#ax = plt.subplot(111)
##ax = midiragnplot(ax, inobssamplefile, survey, save)
#ax = s06_bpt(s06inputfile, s06outputfile, eco, resolve, full, sdsscat, save, ax)

#def totalgals(df, wavelength):
#    return len(df)
#
#def agn(df, wavelength):
#    if wavelength == 'optical':
#        return np.sum(~df['defstarform'])
#    if wavelength == 'midir':
#        return np.sum(df.agnflag)
#
#def totaldwarfs(df, wavelength):
#    bary = pd.read_csv('ECO+RESOLVE_barysample.csv')
#    bary.index = bary.name
#    return np.sum(bary.logmstar.loc[df.name] < 9.5)
#
#def dwarfagn(df, wavelength):
#    if wavelength == 'optical':
#        bary = pd.read_csv('ECO+RESOLVE_barysample.csv')
#        bary.index = bary.name
#        return np.sum(bary.logmstar.loc[df.galname[~df['defstarform']]] < 9.5)
#
#    if wavelength == 'midir':
#        return np.sum(df.logmstar.loc[df.name[df['agnflag']]] < 9.5)
#
#def filenames(index, survey, colname, cat = 'none'):
#    if survey == 'RESOLVE' or survey =='ECO':
#        if 'agn' in colname:
#            files = {'trad': survey.lower()+'_emlineclass_bpt1snr5_'+cat+'.csv', 
#                's06': survey.lower()+'_s06emlineclass_dext_hasnr5_'+cat+'.csv', 
#                'compre': survey.lower()+'_emlineclass_dext_snr5_'+cat+'.csv',
#                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
#                'baryparent': survey+'_barysample.csv'}
#        else:
#            
#            files = {'trad': survey+'_full_bpt1snr5_dext_'+cat+'.csv', 
#                's06': survey+'_full_bpt1snr5_dext_'+cat+'.csv', 
#                'compre': survey+'_full_snr5_dext_'+cat+'.csv',
#                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
#                'baryparent': survey+'_barysample.csv'}
#        
#    else:
#        if 'agn' in colname:
#            files = {'trad': survey.lower()+'_emlineclass_bpt1snr5_'+cat+'.csv', 
#                's06': survey.lower()+'_s06emlineclass_dext_hasnr5_'+cat+'.csv', 
#                'compre': survey.lower()+'_emlineclass_dext_snr5_'+cat+'.csv',
#                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
#                'baryparent': survey+'_barysample.csv'}
#        else:
#            files = {'trad': survey+'_bpt1snr5_dext_'+cat+'.csv', 
#                's06': survey+'_bpt1snr5_dext_'+cat+'.csv', 
#                'compre': survey+'_snr5_dext_'+cat+'.csv',
#                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
#                'baryparent': survey+'_barysample.csv'}
#    
#    return files[index]
#
##printing stats
##surveys = ['ECO+RESOLVE'] #['RESOLVE', 'ECO', 'ECO+RESOLVE']
#surveys = ['RESOLVE']
#sdsscats = ['jhu', 'port', 'nsa']
#columns = ['totalgals', 'agn', 'totaldwarfs', 'dwarfagn']
#waves = ['optical', 'midir']
#
#stats = pd.DataFrame(index = ['trad', 's06', 'compre', 'midir'])#,'baryparent'])
#
#
#for ndx in list(stats.index):
#    for cols in columns:
#        for survey in surveys:
#            if (ndx == 'midir') | (ndx == 'baryparent'):
#                colname = survey+' '+cols
#                df = pd.read_csv(filenames(ndx, survey, colname))
#                if 'name' in df.keys():
#                    df.index = df.name
#                if 'galname' in df.keys():
#                    df.index = df.galname
#
#                wave = 'midir'
#                if colname not in list(stats.keys()):
#                    stats[colname] = np.zeros(len(stats)) 
#                stats[colname].loc[ndx] = eval(cols+'(df,wave)')
#            else:
#                for cat in sdsscats:
#                    colname = survey+' '+cat+' '+cols                    
#                    df = pd.read_csv(filenames(ndx, survey, colname, cat = cat))
#                    if 'name' in df.keys():
#                        df.index = df.name
#                    if 'galname' in df.keys():
#                        df.index = df.galname
#                    wave = 'optical'
#
#                    if colname not in list(stats.keys()):
#                        stats[colname] = np.zeros(len(stats)) 
#                    stats[colname].loc[ndx] = eval(cols+'(df,wave)')
##    
##    
##
##

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
#from bpt_s06 import s06_bpt
#from bpt_plots_new import bpt_plots
#import matplotlib.gridspec as gridspec
#os.chdir("../SDSS_spectra/mid_ir")
#from midir_agn import midiragnplot
#pd.set_option('display.max_columns', None)
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
#he2_flag = 0
#save = 0
#resolve = 0
#eco = 0
#full = 1
##sdsscat = 'port'
#
#sdsscat = 'jhu'
## = 'nsa'
##sdsscat = 'master'
#if sys.platform == 'linux2':
#        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')
#
#else:
#    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
#
#
#if full: 
#    s06inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_bpt1snr5_dext_'+sdsscat+'.csv'#'ECO_full_snr5.csv'
#    bptinputfile = 'ECO+RESOLVE_snr5_dext_'+sdsscat+'.csv' #'ECO_full_bary_jhu.csv'#'ECO_full_snr5.csv'
#    print 'ECO+RESOLVE RESULTS'
#    s06outputfile = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#    bptoutputfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
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
#
##plt.figure()
##gs0 = gridspec.GridSpec(2, 1)
##gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0])
##gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1])
##
##ax1 = plt.subplot(gs00[0, 0])
##ax2 = plt.subplot(gs00[0, 1])
##
##ax3 = plt.subplot(gs01[0, 0])
##ax4 = plt.subplot(gs01[0, 1])
##ax5 = plt.subplot(gs01[0, 2])
#ax1 = plt.subplot(111)
#
#ax1 = s06_bpt(s06inputfile, s06outputfile, eco, resolve, full, sdsscat, save, ax1)
##ax2 = midiragnplot(ax2,inobssamplefile, survey, save)
##os.chdir("C:\Users\mugdhapolimera\github\SDSS_Spectra\/")
##ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile, eco, resolve, full, sdsscat, 
##                          save, ax3, ax4, ax5)
##
##plt.show()
#
##inobssamplefile = bptinputfile
##plt.figure()    
##fig, (ax3, ax4, ax5) = plt.subplots(1,3, sharey = True)
##ax3, ax4, ax5 = bpt_plots(bptinputfile, bptoutputfile, eco, resolve, full, sdsscat, 
##                          save, ax3, ax4, ax5)
##ax = plt.subplot(111)
###ax = midiragnplot(ax, inobssamplefile, survey, save)
##ax = s06_bpt(s06inputfile, s06outputfile, eco, resolve, full, sdsscat, save, ax)
#
##def totalgals(df, wavelength):
##    return len(df)
##
##def agn(df, wavelength):
##    if wavelength == 'optical':
##        return np.sum(~df['defstarform'])
##    if wavelength == 'midir':
##        return np.sum(df.agnflag)
##
##def totaldwarfs(df, wavelength):
##    bary = pd.read_csv('ECO+RESOLVE_barysample.csv')
##    bary.index = bary.name
##    return np.sum(bary.logmstar.loc[df.name] < 9.5)
##
##def dwarfagn(df, wavelength):
##    if wavelength == 'optical':
##        bary = pd.read_csv('ECO+RESOLVE_barysample.csv')
##        bary.index = bary.name
##        return np.sum(bary.logmstar.loc[df.galname[~df['defstarform']]] < 9.5)
##
##    if wavelength == 'midir':
##        return np.sum(df.logmstar.loc[df.name[df['agnflag']]] < 9.5)
##
##def filenames(index, survey, colname, cat = 'none'):
##    if survey == 'RESOLVE' or survey =='ECO':
##        if 'agn' in colname:
##            files = {'trad': survey.lower()+'_emlineclass_bpt1snr5_'+cat+'.csv', 
##                's06': survey.lower()+'_s06emlineclass_dext_hasnr5_'+cat+'.csv', 
##                'compre': survey.lower()+'_emlineclass_dext_snr5_'+cat+'.csv',
##                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
##                'baryparent': survey+'_barysample.csv'}
##        else:
##            
##            files = {'trad': survey+'_full_bpt1snr5_dext_'+cat+'.csv', 
##                's06': survey+'_full_bpt1snr5_dext_'+cat+'.csv', 
##                'compre': survey+'_full_snr5_dext_'+cat+'.csv',
##                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
##                'baryparent': survey+'_barysample.csv'}
##        
##    else:
##        if 'agn' in colname:
##            files = {'trad': survey.lower()+'_emlineclass_bpt1snr5_'+cat+'.csv', 
##                's06': survey.lower()+'_s06emlineclass_dext_hasnr5_'+cat+'.csv', 
##                'compre': survey.lower()+'_emlineclass_dext_snr5_'+cat+'.csv',
##                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
##                'baryparent': survey+'_barysample.csv'}
##        else:
##            files = {'trad': survey+'_bpt1snr5_dext_'+cat+'.csv', 
##                's06': survey+'_bpt1snr5_dext_'+cat+'.csv', 
##                'compre': survey+'_snr5_dext_'+cat+'.csv',
##                'midir': 'mid_ir/'+survey+'_WISE_good.csv',
##                'baryparent': survey+'_barysample.csv'}
##    
##    return files[index]
##
###printing stats
###surveys = ['ECO+RESOLVE'] #['RESOLVE', 'ECO', 'ECO+RESOLVE']
##surveys = ['RESOLVE']
##sdsscats = ['jhu', 'port', 'nsa']
##columns = ['totalgals', 'agn', 'totaldwarfs', 'dwarfagn']
##waves = ['optical', 'midir']
##
##stats = pd.DataFrame(index = ['trad', 's06', 'compre', 'midir'])#,'baryparent'])
##
##
##for ndx in list(stats.index):
##    for cols in columns:
##        for survey in surveys:
##            if (ndx == 'midir') | (ndx == 'baryparent'):
##                colname = survey+' '+cols
##                df = pd.read_csv(filenames(ndx, survey, colname))
##                if 'name' in df.keys():
##                    df.index = df.name
##                if 'galname' in df.keys():
##                    df.index = df.galname
##
##                wave = 'midir'
##                if colname not in list(stats.keys()):
##                    stats[colname] = np.zeros(len(stats)) 
##                stats[colname].loc[ndx] = eval(cols+'(df,wave)')
##            else:
##                for cat in sdsscats:
##                    colname = survey+' '+cat+' '+cols                    
##                    df = pd.read_csv(filenames(ndx, survey, colname, cat = cat))
##                    if 'name' in df.keys():
##                        df.index = df.name
##                    if 'galname' in df.keys():
##                        df.index = df.galname
##                    wave = 'optical'
##
##                    if colname not in list(stats.keys()):
##                        stats[colname] = np.zeros(len(stats)) 
##                    stats[colname].loc[ndx] = eval(cols+'(df,wave)')
###    
###    
###
###
