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

print('SKIF Extended Example for 1-1 # 2:')
print('Extracts electric fields from the files. Plots intensity distributions and spectra.')

#harmonics number
harm1 = 11
harm2 = 13
harm3 = 17
harm4 = 23

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

distance = 20.
a = 10
#**********************Input Parameters:
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
spec1FileName = 'wfr_spec1_1_4.wfr' #for spec1
spec2FileName = 'wfr_spec2_1_4.wfr' #for spec2
wfr2FileName = 'wfr_harm1.wfr' #for harm2
wfr3FileName = 'wfr_harm2.wfr' #for harm3
wfr4FileName = 'wfr_harm3.wfr' #for harm4
wfr5FileName = 'wfr_harm4.wfr' #for harm5
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

x, y= skf.calc_FWHM(wfr1, units='urad')
rms_x, rms_y = skf.calc_bandwidth(wfr1, units='urad')
#print(x,y)
#print(2.335*rms_x, 2.335*rms_y)
#%%
#%%
            ######### Spectrum #######
#skf.skf_plot_spec(spec1)
#skf.skf_plot_spec(spec2)

#%%
            ######### Intensity #######
filepath='/home/andrei/Documents/SKIF_XAS_beamline/1_1/TechReports/inter1/pic'
skf.skf_wfr_subplot_XY(wfr1, fourth_plot=0, save_fig=False, file_path=filepath, figure_name='11_harm_befoure_optics.pdf')
#skf.skf_wfr_subplot_XY(wfr2, fourth_plot=0, save_fig=True, file_path=filepath, figure_name='13_harm_befoure_optics.pdf', show=False)
#skf.skf_wfr_subplot_XY(wfr3, fourth_plot=0, save_fig=True, file_path=filepath, figure_name='17_harm_befoure_optics.pdf', show=False)
#skf.skf_wfr_subplot_XY(wfr4, fourth_plot=0, save_fig=True, file_path=filepath, figure_name='23_harm_befoure_optics.pdf', show=False)
#%% 
    
######### Power density ###########
            
#skf.skf_power_subplot_XY(stkP, units='mm')
#plt.savefig('/home/andrei/Documents/9_term/diplom/beamlines/1_1/power.png', dpi=350)#, bbox_inches='tight')
#

#%%












