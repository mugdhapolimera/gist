# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 13:41:01 2020

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("UVES/gident_437.dat",delimiter= '\s+', \
                 names = ['Sequence','CENTER','INT','FWHM','FLUX'])
plt.figure()
plt.plot(np.array(df['CENTER']),np.array(df['FLUX']))
