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

distance = 25.
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

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr3FileName, 'rb')
wfr3 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr4FileName, 'rb')
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
            ######### Spectrum #######
#skf.skf_plot_spec(spec1)
#skf.skf_plot_spec(spec2)
#
##%%
#save=True
#            ######### Intensity #######
#skf.skf_wfr_subplot_XY(wfr1, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h1+'_harm_befoure_crystal.pdf')
#skf.skf_wfr_subplot_XY(wfr2, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h2+'_harm_befoure_crystal.pdf', show=False)
#skf.skf_wfr_subplot_XY(wfr3, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h3+'_harm_befoure_crystal.pdf', show=False)
#skf.skf_wfr_subplot_XY(wfr4, fourth_plot=0, save_fig=save, file_path=SKIF_path + FigPath, figure_name=h4+'_harm_befoure_crystal.pdf', show=False)
#%% 
  
######### Power density ###########
skf.skf_power_subplot_XY(stkP, wfr=spec1, units='mm')
plt.tight_layout()
plt.savefig(SKIF_path + FigPath + 'power_dens_1-2.pdf')#, bbox_inches='tight')
#%%
plt.figure(figsize=(1.5*4,1.5*4))
wfr=spec1
xy_unit=r'$[mm]$'
fontsize=16
A = 1e3
z = []
Z = []
i = 0 
j = 0 
cmap_ph = plt.get_cmap('viridis')
    
x = np.linspace(A*stkP.mesh.xStart, A*stkP.mesh.xFin, stkP.mesh.nx)
y = np.linspace(A*stkP.mesh.yStart, A*stkP.mesh.yFin, stkP.mesh.ny)
 
for j in range(len(y)):
    for i in range(len(x)):
        z.extend([stkP.arS[j*len(x) + i]])
    Z.append(z)
    z = []
Z = np.array(Z)
plt.pcolormesh(x, y, Z, cmap=cmap_ph)  
plt.ylabel(xy_unit, fontsize=fontsize, labelpad = 0.0, rotation=90)
#plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=18, labelpad = 0.0)
plt.title(r'Плотность мощности', fontsize=fontsize)
plt.ylim(A*stkP.mesh.yStart, A*stkP.mesh.yFin)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xlabel(r'$[мм]$', fontsize=fontsize, labelpad = 0.0)
plt.ylabel(r'$[мм]$', fontsize=fontsize, labelpad = 0.0, rotation=90)
plt.tight_layout()
plt.show()
plt.savefig(SKIF_path + FigPath + 'power_dens_1-2_2d.pdf')#, bbox_inches='tight')

plt.figure(figsize=(1.5*4,1.5*4))
E, spec = skf.renorm_wfr(wfr, emittance=0, elec_fld_units='W/mm^2/eV')
plt.plot(E/1000, spec, color='blue')
plt.xlim(0)
plt.ylim(0)
plt.grid(True)
plt.xlabel(r'$E, [кэВ]$', fontsize=fontsize, labelpad = 0.0)
#    plt.ylabel(r'$I, [\gamma/с/мм^2/0.1\%bw]$', fontsize=fontsize, labelpad = 0.0, rotation=90)
plt.ylabel(r'$I, [Вт/мм^2/эВ]$', fontsize=fontsize, labelpad = 0.0, rotation=90)
#    skf.skf_plot(E, spec, elec_fld_units='ph/s/mm^2/0.1%bw', scale_x = 1000, show=False)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize) 
plt.tight_layout()
plt.show()
plt.savefig(SKIF_path + FigPath + 'Spec_1-2.pdf')#, bbox_inches='tight')

#%%












