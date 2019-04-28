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

print('SKIF Extended Example for 1-2 # 2:')
print('Extracts electric fields from the files. Plots intensity distributions and spectra.')
station = '1_2'

#harmonics number
harm1 = 5
h1 = '5'
harm2 = 7
h2 = '7'
harm3 = 9
h3 = '9'
harm4 = 13
h4= '13'

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

distance = 20.
a = 10
#**********************Input Parameters:
#**********************Input Parameters:
SKIF_path = skf.get_SKIF_directory() #get SKIF project root directory
TablesPath = skf.path_in_project('/' + station + '/TechReports/tabl/')#, your_sys='Mac OC')
FigPath = skf.path_in_project('/' + station + '/TechReports/pic/')
wfrPath = skf.path_in_project('/' + station + '/fields_' + station + '/')
Diamond_T_path = skf.path_in_project('/' + station + '/crystals_data_' + station + '/diamond_T/')

wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
spec1FileName = 'wfr_spec1_' + station + '.wfr' #for spec1
spec2FileName = 'wfr_spec2_' + station + '.wfr' #for spec2
wfr1FileName = 'wfr_harm1.wfr' #for harm2
wfr2FileName = 'wfr_harm2.wfr' #for harm3
wfr3FileName = 'wfr_harm3.wfr' #for harm4
wfr4FileName = 'wfr_harm4.wfr' #for harm5
stkPFileName = 'stkP.wfr'#for power density


#x = np.linspace(-5,5, 10000)
#sigma = 1
#arIx = (1./np.sqrt(2*np.pi*sigma**2))*np.exp(-x**2 / (2*sigma**2))
#plt.plot(x,arIx)
#
#difference = np.max(arIx) - np.min(arIx)
#HM = difference / 2
#nearest = (np.abs(arIx - HM)).argmin()
#print(x[nearest])
#max_index = np.max(np.where(arIx == np.amax(arIx)))#arIx.index(np.max(arIx))
#print(max_index)
#FWHM_x = 2*(np.abs(x[nearest] - x[max_index]))
#print(FWHM_x/2.355)
#%%
print('extracting wfr from files')
afile = open(wfrPathName + spec1FileName, 'rb')
spec1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + spec2FileName, 'rb')
spec2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr3FileName, 'rb')
wfr2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr4FileName, 'rb')
wfr3 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr5FileName, 'rb')
wfr4 =  pickle.load(afile)
afile.close()

#
afile = open(wfrPathName + stkPFileName, 'rb')
stkP =  pickle.load(afile)
afile.close()
print('finishing extracting')

x, y = skf.calc_FWHM(wfr1, units='urad')
rms_x, rms_y = skf.calc_bandwidth(wfr1, units='urad')
#print(x,y)
#print(2.335*rms_x, 2.335*rms_y)
#%%
#%%
            ######### Spectrum #######
#skf.skf_plot_spec(spec1)
#skf.skf_plot_spec(spec2)

#%%
save=True
            ######### Intensity #######
skf.skf_wfr_subplot_XY(wfr1, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h1+'_harm_befoure_crystal.pdf')
skf.skf_wfr_subplot_XY(wfr2, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h2+'_harm_befoure_crystal.pdf', show=False)
skf.skf_wfr_subplot_XY(wfr3, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h3+'_harm_befoure_crystal.pdf', show=False)
skf.skf_wfr_subplot_XY(wfr4, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h4+'_harm_befoure_crystal.pdf', show=False)
#%% 
  
######### Power density ###########
skf.skf_power_subplot_XY(stkP, wfr=spec1, units='urad', save_fig=True, path_name=SKIF_path + FigPath, figure_name='power_dens.pdf')
#%%












