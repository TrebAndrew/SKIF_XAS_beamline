#############################################################################
#Extracts two electric fields from files. Plots intensity distributions and spectra before the filter and after it.
#############################################################################
'''@author Andrei Trebushinin'''

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import SKIF_lib as skf

print('SKIF Extended Example for 1-4:')
print('Extracts electric fields from the files. Plots intensity distributions and spectra.')
station = '1_4'

#harmonics number
harm1 = 3

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

distance = 25.
a = 10
#**********************Input Parameters:
SKIF_path = skf.get_SKIF_directory() #get SKIF project root directory
TablesPath = skf.path_in_project('/' + station + '/TechReports/tabl/')#, your_sys='Mac OC')
FigPath = skf.path_in_project('/' + station + '/TechReports/pic/')
wfrPath = skf.path_in_project('/' + station + '/fields_' + station + '/')
Diamond_T_path = skf.path_in_project('/' + station + '/crystals_data_' + station + '/diamond_T/')

wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
spec1FileName = 'wfr_spec1_' + station + '.wfr' #for spec1
spec2FileName = 'wfr_spec2_' + station + '.wfr' #for spec2
wfr1FileName = 'wfr_harm1.wfr' 
stkPFileName = '25stkP.wfr'#for power density

#%%
print('extracting wfr from files')
afile = open(wfrPathName + spec1FileName, 'rb')
spec1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + spec2FileName, 'rb')
spec2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

#
afile = open(wfrPathName + stkPFileName, 'rb')
stkP =  pickle.load(afile)
afile.close()
print('finishing extracting')

x, y = skf.calc_FWHM(wfr1, units='urad')
rms_x, rms_y = skf.calc_bandwidth(wfr1, units='urad')

#%%
            ######### Spectrum #######
skf.skf_plot_spec(spec1)
skf.skf_plot_spec(spec2)

#%%
            ######### Intensity #######
filepath= SKIF_path + FigPath
skf.skf_wfr_subplot_XY(wfr1, save_fig=True, figure_name='3_harm_before_optics.pdf', fourth_plot=5, show=False, file_path=filepath) #function to draw xy distribution
#%%     
        ######### Power density ###########           
skf.skf_power_subplot_XY(stkP, wfr=spec1, units='mm')
plt.savefig(SKIF_path + FigPath + 'power_dens_1-4.pdf')#, bbox_inches='tight')





