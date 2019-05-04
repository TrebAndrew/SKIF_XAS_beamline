#############################################################################
#SKIF Extended Example for 1-2 # 3:
#Extract electric fields from files. Plot intensity distributions and spectra. 
#Propagates electric field files through CCMs and the DCM
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

print('SKIF Extended Example for 1-2 # 3:')
print('Extracts two electric fields from files. Plots intensity distributions and spectra.\n' 
      'Propagates electric field files through CCMs and the DCM')
station = '1_2'

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
h = 4,135667662e-15
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

#harmonics number
harm1 = 5
h1 = '5'
harm2 = 7
h2 = '7'
harm3 = 9
h3 = '9'
harm4 = 13
h4= '13'
harm = [harm1, harm2, harm3, harm4]

#undulator parameters
Length = 2 # [m]
undper = 0.0156 # [m]
numper = 128
magf = 1.06 # [T]

#**********************Input Parameters:
SKIF_path = skf.get_SKIF_directory() #get SKIF project root directory
TablesPath = skf.path_in_project('/' + station + '/TechReports/tabl/')#, your_sys='Mac OC')
FigPath = skf.path_in_project('/' + station + '/TechReports/pic/')
wfrPath = skf.path_in_project('/' + station + '/fields_' + station + '/')
Diamond_T_path = skf.path_in_project('/' + station + '/crystals_data_' + station + '/diamond_T/')

wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
spec1FileName = 'wfr_spec1_' + station + '.wfr' #for spec1
spec2FileName = 'wfr_spec2_' + station + '.wfr' #for spec2
wfr1FileName = 'wfr_harm1.wfr' #for harm1
wfr2FileName = 'wfr_harm2.wfr' #for harm2
wfr3FileName = 'wfr_harm3.wfr' #for harm3
wfr4FileName = 'wfr_harm4.wfr' #for harm4
stkPFileName = 'stkP.wfr'#for power density

K = 0.9336 * magf * undper * 100 #undulator parameter
E1 = round(4*np.pi*speed_of_light*h_bar*gamma**2/(undper*(1 + K**2/2)), 2) #energy of the first harmonic
for i in range(1, 25, 2):
    Delta_theta = np.sqrt(4*np.pi*speed_of_light*h_bar/(i*E1)/undper/numper) #angle divergence (the first minimum of intensity)
    print('E{} = '.format(i), round(i*E1, 2), '  ang_{} = '.format(i), round(Delta_theta,7))

#%% extracting electric field files from directory wfr PathName
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

absorb = False
#%% Drawing part of the code (before propogation)
print('finishing extracting')
A = 1e6/wfr1.mesh.zStart #scaling facor for the transferring xy distribution from [m] to [rad] 

save=False
filepath= SKIF_path + FigPath
skf.skf_wfr_subplot_XY(wfr1, save_fig=save, figure_name=h1 +'_harm_before_crystal.pdf', fourth_plot=5, show=False, file_path=filepath) #function to draw xy distribution
skf.skf_wfr_subplot_XY(wfr2, save_fig=save, figure_name=h2 +'_harm_before_crystal.pdf', fourth_plot=1, show=False, file_path=filepath) #function to draw xy distribution
skf.skf_wfr_subplot_XY(wfr3, save_fig=save, figure_name=h3 +'_harm_before_crystal.pdf', fourth_plot=1, show=False, file_path=filepath) #function to draw xy distribution
skf.skf_wfr_subplot_XY(wfr4, save_fig=save, figure_name=h4 +'_harm_before_crystal.pdf', fourth_plot=1, show=False, file_path=filepath) #function to draw xy distribution

#%%
wfr = [wfr1, wfr2, wfr3, wfr4]
HARM = []
for (fld, n) in zip(wfr, harm): 
    sigma_x_mm, sigma_y_mm   = skf.calc_bandwidth(fld, units='mm')
    sigma_x_rad, sigma_y_rad = skf.calc_bandwidth(fld, units='urad')
    HARM.append([int(n), sigma_x_mm, sigma_y_mm, sigma_x_rad, sigma_y_rad])
np.savetxt(SKIF_path + TablesPath + "RMS_before.csv", HARM, fmt='%10.d,%10.3f,%10.3f,%10.3f,%10.3f', delimiter=',')#, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')
np.savetxt(SKIF_path + TablesPath + "RMS_before.csv", HARM, fmt='%10.d,%10.3f,%10.3f,%10.3f,%10.3f', delimiter=',')#, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')
#%%
sigma_x, sigma_y = skf.calc_bandwidth(wfr1, units='urad') #calculates sigma fot the electric field distribution 
print('sigma_x_11 = ', round(sigma_x,3),'[urad] \t','sigma_y_11 = ', round(sigma_y,3),'[urad]')
sigma_x, sigma_y = skf.calc_bandwidth(wfr2, units='urad') #calculates sigma fot the electric field distribution 
print('sigma_x_13 = ', round(sigma_x,3),'[urad] \t','sigma_y_13 = ', round(sigma_y,3),'[urad]')
sigma_x, sigma_y = skf.calc_bandwidth(wfr3, units='urad') #calculates sigma fot the electric field distribution 
print('sigma_x_17 = ', round(sigma_x,3),'[urad] \t','sigma_y_17 = ', round(sigma_y,3),'[urad]')
sigma_x, sigma_y = skf.calc_bandwidth(wfr4, units='urad') #calculates sigma fot the electric field distribution 
print('sigma_x_23 = ', round(sigma_x,3),'[urad] \t','sigma_y_23 = ', round(sigma_y,3),'[urad]')

#%%
#skf.skf_plot_spec(spec1) #function to draw spectrum
#skf.skf_plot_spec(spec2) #function to draw spectrum
#%%
#path_name = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/TechReports/inter1/pic/'
#figure_name = 'power_dens.pdf'
#skf.skf_power_subplot_XY(stkP, spec2,units='urad', figure_name=figure_name, path_name=path_name, save_fig=True) #function to power density distribution

#%% A Long long part of the code with the crystal definition

#Diamond(111) Crystal Constants:
dSpC400 = 2.0592929401455575#(Diamond(111))#0.89178 #Crystal reflecting planes d-spacing for Diamond(111) crystal
#psi0rC400 = -0.17530e-04; psi0iC400 = 0.21089e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
psi0rC400 = -6.9764e-6; psi0iC400 = 2.6068e-9#-0.21732e-04; psi0iC400 = 0.28005e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
#psihrC400 = -0.45300e-05; psihiC400 = 0.20314E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihrC400 = 2.5369e-6;   psihiC400 = 1.8151e-9#-0.54377e-05; psihiC400 = 0.25934E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihbrC400 = psihrC400; psihbiC400 = psihiC400 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
thickCryst = 100e-06 #0.5e-03 #Thickness of each crystal [m]
angAsCryst = 0 #Asymmetry angle of each crystal [rad]

#1st Crystal: (for harm 11)

distance = wfr1.mesh.zStart 
opCr1 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)

#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr1 = opCr1.find_orient(wfr1.avgPhotEn)
orientCr1 = orientDataCr1[0] #1st crystal orientation
tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
opCr1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
angle1 = round(np.arctan(nCr1[1]/nCr1[2])*180/np.pi,8)

#%%
#2st Crystal: (for harm 13)

opCr2 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr2 = opCr2.find_orient(wfr3.avgPhotEn)
orientCr2 = orientDataCr2[0] #2st crystal orientation
tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
#Set crystal orientation:
opCr2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])
orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
TrM = [rxCr2, ryCr2, rzCr2] #Input/Output beam transformation matrix (for debugging)
angle2 = round(np.arctan(nCr2[1]/nCr2[2])*180/np.pi,8)
#%%
#3st Crystal: (for harm 23)

opCr3 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr3 = opCr3.find_orient(wfr4.avgPhotEn)
orientCr3 = orientDataCr3[0] #1st crystal orientation
tCr3 = orientCr3[0]; nCr3 = orientCr3[2] # Tangential and Normal vectors to crystal surface
#print('   3st crystal orientation:');# print(' t=', tCr4,'\n', 's=', orientCr4[1],'\n', 'n=', nCr4,'\n')
#Set crystal orientation:
opCr3.set_orient(nCr3[0], nCr3[1], nCr3[2], tCr3[0], tCr3[1])
orientOutFrCr3 = orientDataCr3[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr3 = orientOutFrCr3[0]; ryCr3 = orientOutFrCr3[1]; rzCr3 = orientOutFrCr3[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
TrM = [rxCr3, ryCr3, rzCr3] #Input/Output beam transformation matrix (for debugging)
angle3 = round(np.arctan(nCr3[1]/nCr3[2])*180/np.pi,8)


#%% Double Crystal Monochromator (for harm 17)
###
#C(400) Crystal Constants:
dSpC400 = 3.1355713563754857 #Crystal reflecting planes d-spacing for C(400) crystal
#psi0rC400 = -0.17530e-04; psi0iC400 = 0.21089e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
psi0rC400 = -1.207842005e-5; psi0iC400 = 2.263482755e-7 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
#psihrC400 = -0.45300e-05; psihiC400 = 0.20314E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihrC400 = -6.385703371e-6; psihiC400 = 1.580304013e-7 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihbrC400 = psihrC400; psihbiC400 = psihiC400 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
thickCryst = 10.0e-03 #0.5e-03 #Thickness of each crystal [m]
angAsCryst = 0 #Asymmetry angle of each crystal [rad]

#1st Crystal:
DCM_Cr1 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr1 = DCM_Cr1.find_orient(wfr2.avgPhotEn)
orientCr1 = orientDataCr1[0] #1st crystal orientation
tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
#print('   1st crystal orientation:'); print('   t=', tCr1, 's=', orientCr1[1], 'n=', nCr1)
#Set crystal orientation:
DCM_Cr1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
angleDCM_Cr1 = round(np.arctan(nCr1[1]/nCr1[2])*180/np.pi,8)
orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
uti_math.matr_print(TrM)

#2nd Crystal:
DCM_Cr2 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst)
#Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr2 = DCM_Cr2.find_orient(wfr2.avgPhotEn, _ang_dif_pl=3.141593)
orientCr2 = orientDataCr2[0] #2nd crystal orientation
tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
#print('   2nd crystal orientation:'); print('   t=', tCr2, 's=', orientCr2[1], 'n=', nCr2)
#Set crystal orientation:
DCM_Cr2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])
orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   2nd crystal output beam frame:'); print('   ex=', rxCr2, 'ey=', ryCr2, 'ez=', rzCr2)
TrM = uti_math.matr_prod(TrM, [rxCr2, ryCr2, rzCr2]) #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
uti_math.matr_print(TrM)
#print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')


#%%
#***************** Optical Elements and Propagation Parameters
distance = wfr1.mesh.zStart 
sigma_x_mm, sigma_y_mm   = skf.calc_bandwidth(wfr1, units='mm')
slit_x = 2*sigma_x_mm#4*sigma_x*distance*1e-3 # mm
slit_y = 2*sigma_y_mm

#Drift_AFTER_Cr   = SRWLOptD(1.) #Drift from diamond filter
Drift_BEFORE_Cr1 = SRWLOptD(1.) #Drift from the first diamond(111) CCM for 5 harmonic
Drift_BEFORE_Cr2 = SRWLOptD(7.) #Drift before DCM for 7 harmonic
Drift_BEFORE_Cr3 = SRWLOptD(2.) #Drift before the second diamond(111) CCM for 9 harmonic
Drift_BEFORE_Cr4 = SRWLOptD(4.) #Drift before the third diamond(111) CCM for 13 harmonic

SSA = SRWLOptA('r', 'a', 2*slit_x*1e-03, 2*slit_y*1e-03) #slit parameter
Drift_AFTER_SLIT = SRWLOptD(1.) #Drift after the slit
#Drift_SSA_SCREEN = SRWLOptD(1.) #Drift from SSA to Center of Vertically Focusing K-B Mirror (VKB)

#Wavefront Propagation Parameters:
#                       [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
prParInit =             [ 0,  0,  1.,  1,  0, 1.0, 2.0, 1.0, 2.0,  0,  0,   0]
pprPar0 =               [ 0,  0,  1.,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppSSA               =   [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppDrift_BEFORE_Cr1  =   [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_BEFORE_Cr2  =   [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_BEFORE_Cr3  =   [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_BEFORE_Cr4  =   [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

prParPost =             [ 0,  0,  1.,  0,  0, 1.0, 3.0, 1.0, 3.0,  0,  0,   0]

ppDrift_AFTER_Cr1  =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_BEFORE_SLIT =   [ 0,  0, 1.0,  1,  0, 1.,  0.7, 1.,   0.7, 0,  0,   0]
ppDrift_AFTER_SLIT  =   [ 0,  0, 1.0,  1,  0, 1.1, 1.0, 1.1,  1.0, 0,  0,   0]

#[ 0]: Auto-Resize (1) or not (0) Before propagation
#[ 1]: Auto-Resize (1) or not (0) After propagation
#[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[ 6]: Horizontal Resolution modification factor at Resizing
#[ 7]: Vertical Range modification factor at Resizing
#[ 8]: Vertical Resolution modification factor at Resizing
#[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
#optBL1 = SRWLOptC([SSA,Drift_BEFORE_Cr1,SSA, Drift_BEFORE_Cr1, opCr1], 
#                  [prParInit, ppSSA,ppDrift_BEFORE_Cr1, ppSSA, ppDrift_BEFORE_Cr1, pprPar0, prParPost])

optBL1 = SRWLOptC([SSA, Drift_BEFORE_Cr1, opCr1], # 5 harm
                  [prParInit, ppSSA, ppDrift_BEFORE_Cr1, pprPar0, prParPost])

optBL2 = SRWLOptC([SSA, Drift_BEFORE_Cr2, DCM_Cr1, DCM_Cr2], # 7 harm
                  [prParInit, ppSSA, ppDrift_BEFORE_Cr2, pprPar0, pprPar0, prParPost])

optBL3 = SRWLOptC([SSA, Drift_BEFORE_Cr3, opCr2], # 9 harm
                  [prParInit, ppSSA, ppDrift_BEFORE_Cr3, pprPar0, prParPost])

optBL4 = SRWLOptC([SSA, Drift_BEFORE_Cr4, opCr3], # 13 harm
                  [prParInit, ppSSA, ppDrift_BEFORE_Cr4, pprPar0, prParPost])
#optBL = SRWLOptC([Drift_BEFORE_SLIT], [ppDrift_BEFORE_SLIT, prParPost])

#%%
#///////////WAVEFRONT PROPOGATION//////#
print('   Simulating Electric Field Wavefront Propagation ... ', end='')
t0 = time.time();
#srwl.PropagElecField(spec1, optBL)
srwl.PropagElecField(wfr1, optBL1)
srwl.PropagElecField(wfr2, optBL2)
srwl.PropagElecField(wfr3, optBL3)
srwl.PropagElecField(wfr4, optBL4)
print('done; lasted', round(time.time() - t0), 's')
#%%
wfrContainer = [wfr1, wfr2, wfr3, wfr4]
angleContainer = [angle1, angleDCM_Cr1, angle2, angle3]
EffCrThick = []
projContainer = []
harmContainer = [harm1, harm2, harm3, harm4]

for (wfr, angle, harm) in zip(wfrContainer, angleContainer, harmContainer):
    sigma_x, sigma_y = skf.calc_bandwidth(wfr, units='mm')
    print(harm, 'harmonic ','$angle = ', round(angle, 2), r'[deg]$', 'at energy $' + str(wfr.mesh.eStart) + '$ eV' + '\\\\', '\n', 
    '$Cr_{effective}  = ', round((100 / np.cos(angle*np.pi/180))), r'[\mu m]$\\',)
    '$proj_x  = ', round((sigma_x / np.cos(angle*np.pi/180)),3), r'[mm]$\\','\n',
    EffCrThick.append(round((100 / np.cos(angle*np.pi/180))))
    projContainer.append(round((sigma_x / np.cos(angle*np.pi/180)),3))

print(EffCrThick)
print(projContainer)
#angleContainer.append(angleDCM_Cr1)
#harmContainer.append(harm3)
EffCrThick[1] = 0
#Save angles to a file
#wfr =  [wfr1, wfr2, wfr3, wfr4]
#harm = [harm1, harm2, harm3, harm4]
HARM = []
    
for (angle, n, thick, proj) in zip(angleContainer, harmContainer, EffCrThick, projContainer): 
    HARM.append([int(n), angle, thick, proj])
np.savetxt(SKIF_path + TablesPath + 'Cr_angles.csv', HARM, fmt='%10.d,%10.3f,%10.d,%10.3f', delimiter=',')#, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')


#%% Drawing part of the code (after propogation)
save=False
filepath = SKIF_path + FigPath
skf.skf_wfr_subplot_XY(wfr1, save_fig=save, figure_name=h1+'_harm_after_crystal.pdf', fourth_plot=5, show=False, file_path=filepath)
skf.skf_wfr_subplot_XY(wfr2, save_fig=save, figure_name=h2+'_harm_after_crystal.pdf', fourth_plot=1, show=False, file_path=filepath)
skf.skf_wfr_subplot_XY(wfr3, save_fig=save, figure_name=h3+'_harm_after_crystal.pdf', fourth_plot=1, show=False, file_path=filepath)
skf.skf_wfr_subplot_XY(wfr4, save_fig=save, figure_name=h4+'_harm_after_crystal.pdf', fourth_plot=1, show=False, file_path=filepath)
#%%
wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
wfr1FileName = 'wfr_harm1_after_Cr.wfr' #for harm2
wfr2FileName = 'wfr_harm2_after_Cr.wfr' #for harm3
wfr3FileName = 'wfr_harm3_after_Cr.wfr' #for harm4
wfr4FileName = 'wfr_harm4_after_Cr.wfr' #for harm5
wfrFileName = [wfr1FileName, wfr2FileName, wfr3FileName, wfr4FileName]#, stkPFileName]
wfrContainer = [wfr1, wfr2, wfr3, wfr4]

print('saving to the files')
#*****************Saving to files
for (wfr, fname) in zip(wfrContainer, wfrFileName):
    print(wfr, fname, "\n")
    afile = open(wfrPathName + fname, 'wb')
    pickle.dump(wfr, afile)
    afile.close()
    
