#############################################################################
#Extract two electric fields from files. Plot intensity distribution and spectrum.
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle
import random as rn
import numpy as np
import matplotlib.pyplot as plt
import SKIF_lib as skf


print('SKIF Extended Example for 1-1 # 3:')
print('Extract two electric fields from files. Plot intensity distribution and spectrum. Propagate a wfr through a CCM')
speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03
e_ = 1.60218e-19 #elementary charge

#harmonics number
harm2 = 11
harm3 = 13
harm4 = 17
harm5 = 23
#undulator parameters
Length = 2.3 # m
undper = 0.018 # m
numper = 128
magf = 1.36  
#**********************Input Parameters:
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
wfr1FileName = 'spec.wfr' # for spectrum
fullspecFileName = 'full_spec.wfr' # for spectrum
wfr2FileName = 'wfr_harm2.wfr' #for harm2
wfr3FileName = 'wfr_harm3.wfr' #for harm3
wfr4FileName = 'wfr_harm4.wfr' #for harm4
wfr5FileName = 'wfr_harm5.wfr' #for harm5

K = 0.9336 * magf * undper * 100
E1 = round(4*np.pi*speed_of_light*h_bar*gamma**2/(undper*(1 + K**2/2)), 2)
for i in range(1, 25, 2):
    Delta_theta = np.sqrt(4*np.pi*speed_of_light*h_bar/(i*E1)/undper/numper)
    print('E{} = '.format(i), round(i*E1, 2), '  ang_{} = '.format(i), round(Delta_theta,7))

#%%
print('extracting wfr from files')
afile = open(wfrPathName + fullspecFileName, 'rb')
wfr0 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr2   =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr3FileName, 'rb')
wfr3 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr4FileName, 'rb')
wfr4 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr5FileName, 'rb')
wfr5 =  pickle.load(afile)
afile.close()
print('finishing extracting')
A = 1e6/wfr2.mesh.zStart

skf.skf_subplot_XY(wfr2, save_fig=True, figure_name='11_harm_befour_crystal')
skf.skf_subplot_XY(wfr3, save_fig=True, figure_name='13_harm_befour_crystal')
skf.skf_subplot_XY(wfr4, save_fig=True, figure_name='17_harm_befour_crystal')
skf.skf_subplot_XY(wfr5, save_fig=True, figure_name='23_harm_befour_crystal')
sigma_x, sigma_y = skf.calc_bandwidth(wfr2)
print('sigma_x_11 = ', round(sigma_x,3),'[urad] \t','sigma_y_11 = ', round(sigma_y,3),'[urad]')

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... \n', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
#uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

#skf.skf_plot_spec(wfr1, crystal=True)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)

E, spec = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV')
skf.skf_plot(E, spec, elec_fld_units='W/mm^2/eV', color='blue')

T = skf.pycry_trans(crystal='diamond', Emin=wfr1.mesh.eStart, Emax=wfr1.mesh.eFin, ne=wfr1.mesh.ne)
spec = spec[:len(T)]
E = E[:len(T)]
#plt.savefig(path_name + figure_name, dpi=150)#, bbox_inches='tight')

srwl.CalcIntFromElecField(arI1, wfr0, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
E2, full_spec = skf.renorm_wfr(wfr0, elec_fld_units='W/mm^2/eV')
skf.skf_plot(E2, full_spec, color='pink', elec_fld_units='W/mm^2/eV', grid=True, linewidth=0.3)

#%%

plt.figure(figsize=(1.5*10,1.5*3))
skf.skf_plot(E, spec*T, color='red', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2)
W = (np.sum(spec) - np.sum(spec*T))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

skf.skf_plot(E, spec*T**2, color='blue', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2)
W1 = ((np.sum(spec*T) - np.sum(spec*T**2)))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

skf.skf_plot(E, spec*T**3, color='red', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2)
W2 = (np.sum(spec*T**2) - np.sum(spec*T**3))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

skf.skf_plot(E, spec*T**4, color='blue', elec_fld_units='W/mm^2/eV', grid=True, linewidth=2)
W3 = (np.sum(spec*T**3) - np.sum(spec*T**4))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

skf.skf_plot(E, spec*T**5, color='black', elec_fld_units='W/mm^2/eV', grid=True, linewidth=3)
W4 = (np.sum(spec*T**4) - np.sum(spec*T**5))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

full_W = np.sum(full_spec)*((wfr0.mesh.eFin - wfr0.mesh.eStart) / wfr0.mesh.ne)
left_W = (np.sum(spec*T**5))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)

t = (r"$W_F = {} \; [Вт/мм^2]$".format(round(full_W,1)) + "\n"
     r"$W_a = {} \; [Вт/мм^2]$".format(round(W,1)) + "\n"
     r"$W1_a = {} \; [Вт/мм^2]$".format(round(W1,1)) + "\n"
     r"$W2_a = {} \; [Вт/мм^2]$".format(round(W2,1)) + "\n"
     r"$W3_a = {} \; [Вт/мм^2]$".format(round(W3,1)) + "\n"
     r"$W4_a = {} \; [Вт/мм^2]$".format(round(W4,1)) + "\n"
     r"$W_L = {} \; [Вт/мм^2]$".format(round(left_W,1)) + "\n")
plt.text(0.7, 0.79, t,
         {'color': 'black', 'fontsize': 14},
         horizontalalignment='left',
         verticalalignment='center',
         rotation=0,
         clip_on=False,
         transform=plt.gca().transAxes)
skf.skf_plot(E2, full_spec, color='pink', elec_fld_units='W/mm^2/eV', grid=True, linewidth=0.3, save_fig=True, file_name='full_spec')

plt.show()

#%%
#Diamond(111) Crystal Constants:
dSpC400 = 2.0592929401455575#(Diamond(111))#0.89178 #Crystal reflecting planes d-spacing for Diamond(111) crystal
#psi0rC400 = -0.17530e-04; psi0iC400 = 0.21089e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
psi0rC400 = -6.9764e-6; psi0iC400 = 2.6068e-9#-0.21732e-04; psi0iC400 = 0.28005e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
#psihrC400 = -0.45300e-05; psihiC400 = 0.20314E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihrC400 = 2.5369e-6;   psihiC400 = 1.8151e-9#-0.54377e-05; psihiC400 = 0.25934E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihbrC400 = psihrC400; psihbiC400 = psihiC400 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
thickCryst = 100e-06 #0.5e-03 #Thickness of each crystal [m]
angAsCryst = 0 #Asymmetry angle of each crystal [rad]
#%%
distance = wfr2.mesh.zStart 
#1st Crystal:
opCr1 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):

orientDataCr1 = opCr1.find_orient(wfr2.avgPhotEn)
orientCr1 = orientDataCr1[0] #1st crystal orientation
tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
print('   1st crystal orientation:');# print(' t=', tCr1,'\n', 's=', orientCr1[1],'\n', 'n=', nCr1,'\n')
#Set crystal orientation:
opCr1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
angle1 = round(np.arctan(nCr1[1]/nCr1[2])*180/np.pi,8)
sigma_x_11,sigma_y_11 = skf.calc_bandwidth(wfr2)
print('   angle_11 = ', angle1, r'[deg]','\n', 
      'sigma_x_11 = ', sigma_x_11*(distance+1)*1e-3, r'[mm]','\n', 
      'beam_proj =', round(abs(sigma_x_11*(distance+1)*1e-3/np.cos(angle1*np.pi/180)),2), r'[mm]','\n')

orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
#uti_math.matr_print(TrM)
#%%

#2st Crystal:
opCr2 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):

orientDataCr2 = opCr2.find_orient(wfr3.avgPhotEn)
orientCr2 = orientDataCr2[0] #1st crystal orientation
tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
print('   2st crystal orientation:');# print(' t=', tCr2,'\n', 's=', orientCr2[1],'\n', 'n=', nCr2,'\n')
#Set crystal orientation:
opCr2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])

angle2 = round(np.arctan(nCr2[1]/nCr2[2])*180/np.pi,8)
sigma_x_13,sigma_y_13 = skf.calc_bandwidth(wfr3)
print('   angle_13 = ', angle2, r'[deg]','\n', 
      'sigma_x_13 = ', sigma_x_13*(distance+2)*1e-3, r'[mm]','\n', 
      'beam_proj =', round(abs(sigma_x_13*(distance+2)*1e-3/np.cos(angle2*np.pi/180)),2), r'[mm]','\n')

orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr2, ryCr2, rzCr2] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
#uti_math.matr_print(TrM)
#%%
#3st Crystal:

opCr3 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):

orientDataCr3 = opCr3.find_orient(wfr4.avgPhotEn)
orientCr3 = orientDataCr3[0] #1st crystal orientation
tCr3 = orientCr3[0]; nCr3 = orientCr3[2] # Tangential and Normal vectors to crystal surface
print('   3st crystal orientation:');# print(' t=', tCr3,'\n', 's=', orientCr3[1],'\n', 'n=', nCr3,'\n')
#Set crystal orientation:
opCr3.set_orient(nCr3[0], nCr3[1], nCr3[2], tCr3[0], tCr3[1])
angle3 = round(np.arctan(nCr3[1]/nCr3[2])*180/np.pi,8)
sigma_x_17,sigma_y_17 = skf.calc_bandwidth(wfr4)
print('   angle_13 = ', angle3, r'[deg]','\n', 
      'sigma_x_17 = ', sigma_x_17*(distance+3)*1e-3, r'[mm]','\n', 
      'beam_proj =', round(abs(sigma_x_17*(distance+3)*1e-3/np.cos(angle3*np.pi/180)),2), r'[mm]','\n')

orientOutFrCr3 = orientDataCr3[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr3 = orientOutFrCr3[0]; ryCr3 = orientOutFrCr3[1]; rzCr3 = orientOutFrCr3[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr3, ryCr3, rzCr3] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
#uti_math.matr_print(TrM)
#%%
#4st Crystal:

opCr4 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=1)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):

orientDataCr4 = opCr4.find_orient(wfr5.avgPhotEn)
orientCr4 = orientDataCr4[0] #1st crystal orientation
tCr4 = orientCr4[0]; nCr4 = orientCr4[2] # Tangential and Normal vectors to crystal surface
print('   4st crystal orientation:');# print(' t=', tCr4,'\n', 's=', orientCr4[1],'\n', 'n=', nCr4,'\n')
#Set crystal orientation:
opCr4.set_orient(nCr4[0], nCr4[1], nCr4[2], tCr4[0], tCr4[1])
angle4 = round(np.arctan(nCr4[1]/nCr4[2])*180/np.pi,8)
sigma_x_23,sigma_y_23 = skf.calc_bandwidth(wfr5)
print('   angle_23 = ', angle4, r'[deg]','\n', 
      'sigma_x_23 = ', sigma_x_17*(distance+4)*1e-3, r'[mm]','\n', 
      'beam_proj =', round(abs(sigma_x_23*(distance+4)*1e-3/np.cos(angle4*np.pi/180)),2), r'[mm]','\n')

orientOutFrCr4 = orientDataCr4[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr4 = orientOutFrCr4[0]; ryCr4 = orientOutFrCr4[1]; rzCr4 = orientOutFrCr4[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr4, ryCr4, rzCr4] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
#uti_math.matr_print(TrM)

#%%
#***************** Optical Elements and Propagation Parameters
distance = wfr2.mesh.zStart 
slit = 4*sigma_x*distance*1e-3 # mm
Drift_AFTER_Cr  = SRWLOptD(1.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)
Drift_BEFORE_Cr1 = SRWLOptD(1.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)
Drift_BEFORE_Cr2 = SRWLOptD(2.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)
Drift_BEFORE_Cr3 = SRWLOptD(3.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)
Drift_BEFORE_Cr4 = SRWLOptD(4.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)

SSA = SRWLOptA('r', 'a', 2*slit*1e-03, 2*slit*1e-03)
Drift_AFTER_SLIT = SRWLOptD(1.)
#Drift_SSA_SCREEN = SRWLOptD(1.) #Drift from SSA to Center of Vertically Focusing K-B Mirror (VKB)

#Wavefront Propagation Parameters:
#                       [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
prParInit =             [ 0,  0,  1.,  1,  0, 1.0, 2.0, 1.0, 2.0,  0,  0,   0]
pprPar0 =               [ 0,  0,  1.,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppSSA               =   [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]

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
#optBL = SRWLOptC([opCr1, Drift_BEFORE_SLIT, SSA, Drift_AFTER_SLIT], [ppDrift_BEFORE_SLIT, ppSSA, ppDrift_AFTER_SLIT, ppopCr1])
optBL1 = SRWLOptC([SSA, Drift_BEFORE_Cr1, opCr1], 
                  [ppSSA, ppDrift_BEFORE_Cr1, pprPar0, prParPost])

optBL2 = SRWLOptC([SSA,Drift_BEFORE_Cr2, opCr2], 
                  [ppSSA, ppDrift_BEFORE_Cr2, pprPar0, prParPost])

optBL3 = SRWLOptC([SSA,Drift_BEFORE_Cr3, opCr3], 
                  [ppSSA, ppDrift_BEFORE_Cr3, pprPar0, prParPost])

optBL4 = SRWLOptC([SSA,Drift_BEFORE_Cr4, opCr4], 
                  [ppSSA, ppDrift_BEFORE_Cr4, pprPar0, prParPost])


#optBL = SRWLOptC([Drift_BEFORE_SLIT], [ppDrift_BEFORE_SLIT, prParPost])

#%%
#///////////WAVEFRONT PROPOGATION//////#
print('   Simulating Electric Field Wavefront Propagation ... ', end='')
t0 = time.time();
#srwl.PropagElecField(wfr1, optBL)
srwl.PropagElecField(wfr2, optBL1)
srwl.PropagElecField(wfr3, optBL2)
srwl.PropagElecField(wfr4, optBL3)
srwl.PropagElecField(wfr5, optBL4)
print('done; lasted', round(time.time() - t0), 's')

skf.skf_subplot_XY(wfr2, save_fig=True, figure_name='11_harm_after_crystal')

skf.skf_subplot_XY(wfr3, save_fig=True, figure_name='13_harm_after_crystal')

skf.skf_subplot_XY(wfr4, save_fig=True, figure_name='17_harm_after_crystal')

skf.skf_subplot_XY(wfr5, save_fig=True, figure_name='23_harm_after_crystal')
#








