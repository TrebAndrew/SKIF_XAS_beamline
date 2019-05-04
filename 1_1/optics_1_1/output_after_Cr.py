#############################################################################
#SKIF Extended Example for 1-1 # 3:
#Extract electric fields from files. Plot intensity distributions and spectra. 
#Propagates electric field files through CCMs(implemented) and the DCM(not implemented yet)
#############################################################################
'''@author Andrei Trebushinin'''

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import SKIF_lib as skf #special labrary for the SKIF project. Add the directory to your PYTHONPATH
import platform


print('SKIF Extended Example for 1-1 # 3:')
print('Extracts two electric fields from files. Plots intensity distributions and spectra.\n' 
      'Propagates electric field files through CCMs and the DCM')
station = '1_1'

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
h = 4.135667662e-15 #[eV*s]
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

#harmonics number
harm1 = 11
harm2 = 13
harm3 = 17
harm4 = 23
harm = [harm1, harm2, harm3, harm4]

#undulator parameters
Length = 2.3 # m
undper = 0.018 # m
numper = 128
magf = 1.36  
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
wfr1FileName = 'wfr_harm1.wfr' #for harm 11 for CCD
wfr2FileName = 'wfr_harm2.wfr' #for harm 13 for CCD
wfr3FileName = 'wfr_harm3.wfr' #for harm 17 for DCM
wfr4FileName = 'wfr_harm4.wfr' #for harm 23 for CCD
stkPFileName = 'stkP.wfr'#for power density

wfrFileName = [wfr1FileName, wfr2FileName, wfr3FileName, wfr4FileName]#, stkPFileName]

K = 0.9336 * magf * undper * 100 #undulator parameter
E1 = round(4*np.pi*speed_of_light*h_bar*gamma**2/(undper*(1 + K**2/2)), 2) #energy of the first harmonic
for i in range(1, 25, 2):
    Delta_theta = np.sqrt(4*np.pi*speed_of_light*h_bar/(i*E1)/undper/numper) #angle divergence (the first minimum of intensity)
    print('E{} = '.format(i), round(i*E1, 2), '  ang_{} = '.format(i), round(Delta_theta,7))

#%%
#extracting electric field files from directory wfr PathName
print('extracting wfr from files')

afile = open(wfrPathName + spec1FileName, 'rb')
spec1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + spec2FileName, 'rb')
spec2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1   =  pickle.load(afile)
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

afile = open(wfrPathName + stkPFileName, 'rb')
stkP =  pickle.load(afile)
afile.close()

wfrContainer = [wfr1, wfr2, wfr3, wfr4]
#%% Drawing spectrum absorption characteristic

plt.figure(figsize=(1.5*10,1.5*3))
z1= spec1.mesh.zStart

E, spec = skf.renorm_wfr(spec1, elec_fld_units='W/mm^2/eV')
skf.skf_plot(E, spec, elec_fld_units='W/mm^2/eV', color='pink', grid=True, linewidth=1.5, show=False)

E2, full_spec = skf.renorm_wfr(spec1, elec_fld_units='W/mm^2/eV')
skf.skf_plot(E2, full_spec, color='pink', elec_fld_units='W/mm^2/eV', grid=True, linewidth=1, show=False)

file_path = SKIF_path + Diamond_T_path
file_name = 'diamond_T_100.txt'
X, T = skf.read_crystal_trans(file_path, file_name)
X_100, T_100 = skf.read_crystal_trans(file_path,file_name='diamond_T_100.txt')
X_480, T_480 = skf.read_crystal_trans(file_path,file_name='diamond_T_480.txt')
X_568, T_568 = skf.read_crystal_trans(file_path,file_name='diamond_T_568.txt')
X_742, T_742 = skf.read_crystal_trans(file_path,file_name='diamond_T_742.txt')
X_1004, T_1004 = skf.read_crystal_trans(file_path,file_name='diamond_T_1004.txt')

#skf.skf_plot(E, spec*T_100, color='red', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2, show=False)
W = (np.sum(spec) - np.sum(spec*T_100))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)

#skf.skf_plot(E, spec*(T_100*T_480)*(z1/(z1+1))**2, color='blue', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2, show=False)
W1 = (np.sum(spec*T_100*(z1/(z1+1))**2) - np.sum(spec*T_100*T_480*(z1/(z1+1))**2))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)

#skf.skf_plot(E, spec*(T_100*T_480*T_568)*(z1/(z1+2))**2, color='red', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2, show=False)
W2 = (np.sum(spec*(T_100*T_480)*(z1/(z1+2))**2) - np.sum(spec*(T_100*T_480*T_568)*(z1/(z1+2))**2))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)

skf.skf_plot(E, spec*(T_100*T_480*T_568*T_742)*(z1/(z1+3))**2, color='black', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2, show=False)
W3 = (np.sum(spec*(T_100*T_480*T_568)*(z1/(z1+3))**2) - np.sum(spec*(T_100*T_480*T_568*T_742)*(z1/(z1+3))**2))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)

#skf.skf_plot(E, spec*(T_100*T_480*T_568*T_742*T_1004)*(z1/(z1+4))**2, color='black', elec_fld_units='W/mm^2/eV', 
#             grid=True, linewidth=3, show=False)
#W4 = (np.sum(spec*T_100*T_480*T_568*T_742*(z1/(z1+4))**2) - np.sum(spec*T_100*T_480*T_568*T_742*T_1004*(z1/(z1+4))**2))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)

full_W = np.sum(full_spec)*((spec2.mesh.eFin - spec2.mesh.eStart) / spec2.mesh.ne)
left_W = (np.sum(spec*T_100*T_480*T_568*T_742*(z1/(z1+4))**2))*((spec1.mesh.eFin - spec1.mesh.eStart) / spec1.mesh.ne)


t = (r"$W_F \approx {} \; [Вт/мм^2]\; падающая\; на\; алмазное\; окно$".format(round(full_W,1)) + "\n"
     r"$W_a \approx {} \; [Вт/мм^2]\; поглощенная\; в\; алмазном\; окне$".format(round(W,1)) + "\n"
     r"$W1_a \approx {} \; [Вт/мм^2] \; поглощенная\; в\; Cr1$".format(round(W1,1)) + "\n"
     r"$W2_a \approx {} \; [Вт/мм^2] \; поглощенная\; в\; Cr2$".format(round(W2,1)) + "\n"
     r"$W3_a \approx {} \; [Вт/мм^2] \; поглощенная\; в\; Cr3$".format(round(W3,1)) + "\n"
     #r"$W4_a \approx {} \; [Вт/мм^2] \; поглощенная\; в\; CCM4$".format(round(W4,1)) + "\n"
     r"$W_L \approx {} \; [Вт/мм^2] \; падающая\; на\; DCM$".format(round(left_W,1)) + "\n") 

plt.text(0.45, 0.59, t,
         {'color': 'black', 'fontsize': 18},
         horizontalalignment='left',
         verticalalignment='center',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
plt.ylim(0, 1.2)
plt.xlim(0, 60000)
plt.savefig(SKIF_path + FigPath + 'spec.pdf', dpi=250)#, bbox_inches='tight')
plt.show()

#%%Parameters extraction

#table of the RMS size of the beams after crystals 
wfr = [wfr1, wfr2, wfr3, wfr4]
harm = [harm1, harm2, harm3, harm4]
HARM = []
for (fld, n) in zip(wfr, harm): 
    sigma_x_mm, sigma_y_mm   = skf.calc_bandwidth(fld, units='mm')
    sigma_x_rad, sigma_y_rad = skf.calc_bandwidth(fld, units='urad')
    HARM.append([int(n), sigma_x_mm, sigma_y_mm, sigma_x_rad, sigma_y_rad])
np.savetxt(SKIF_path + TablesPath + "RMS_after.csv", HARM, fmt='%10.d,%10.3f,%10.3f,%10.3f,%10.3f', delimiter=',')#, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')

#%% Load angles of the Cr
Monchr_ang = np.loadtxt(SKIF_path + TablesPath + "Cr_angles.csv", delimiter=',')#grad
#print(Monchr_ang)
#%%
Darvin_curve_diamond = np.loadtxt(SKIF_path + TablesPath + "Darvin_curve_diamond.csv", delimiter=',')#urad
print(Darvin_curve_diamond)
#%%
Darvin_curve_selicon = np.loadtxt(SKIF_path + TablesPath + "Darvin_curve_silicon.csv", delimiter=',')#urad
print(Darvin_curve_selicon)
#%%Save beam parameter to a file E, lambda, dE/E, RMS, Flux(ph/s), Flux(ph/s/0.1%bw)
wfr = [wfr1, wfr2, wfr3, wfr4]
harm = [harm1, harm2, harm3, harm4]
HARM = []
dE_E=4*[0]
dE_E[0] = Darvin_curve_diamond[0][1]*1e-6/np.abs(np.tan((90.0-Monchr_ang[0][1])*np.pi/180))#relative bandwidth
dE_E[1] = Darvin_curve_diamond[1][1]*1e-6/np.abs(np.tan((90.0-Monchr_ang[1][1])*np.pi/180))#relative bandwidth
dE_E[2] = Darvin_curve_diamond[2][1]*1e-6/np.abs(np.tan((90.0-Monchr_ang[2][1])*np.pi/180))#relative bandwidth
dE_E[3] = Darvin_curve_selicon[0][1]*1e-6/np.abs(np.tan((90.0-Monchr_ang[0][1])*np.pi/180))#relative bandwidth

for (fld, n, bwth) in zip(wfr, harm, dE_E):
    sigma_x_mm, sigma_y_mm   = skf.calc_bandwidth(fld, units='mm')
    sigma_x_rad, sigma_y_rad = skf.calc_bandwidth(fld, units='urad')
#    print(angle[1])
    #dE_E = fwhm[1]*1e-6/np.abs(np.tan((90.0-angle[1])*np.pi/180))#relative bandwidth
    
    E = np.linspace(fld.mesh.eStart, fld.mesh.eFin, fld.mesh.ne) #energy range
    arI = array('f', [0]*fld.mesh.nx*fld.mesh.ny) #"flat" 2D array to take intensity data
    srwl.CalcIntFromElecField(arI, fld, 6, 1, 3, fld.mesh.eStart, 0, 0)
    
    F = np.sum(arI)*(fld.mesh.yFin - fld.mesh.yStart)*(fld.mesh.yFin - fld.mesh.yStart)/fld.mesh.nx/fld.mesh.ny#ph/s/0.1bw
    print('F = ', F)
    arI = arI/E/1e-3 #ph/s/eV
    TotF = np.sum(arI)*bwth*fld.avgPhotEn*(fld.mesh.yFin - fld.mesh.yStart)*(fld.mesh.yFin - fld.mesh.yStart)/fld.mesh.nx/fld.mesh.ny
    print('TotF = ', TotF)
    HARM.append([int(n), fld.avgPhotEn,  h*speed_of_light*1e9/fld.avgPhotEn, TotF, F, bwth])
np.savetxt(SKIF_path + TablesPath +  "ph_beam_par_after_cr.csv", 
           HARM, fmt='%10.d,%10.0f,%10.4f,%.2e,%.2e,%.2e', delimiter=',')#, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')
#%%










