#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:26:32 2019

@author: andrei
"""

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

print('1-4 beamline')
print('Create an undulator for 1-4 station.')
station = 'MechMov'

#**********************Output files
speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03 #relative energy E_electron/m_e [GeV/Gev]
e_ = 1.60218e-19 #elementary charge

#**********************Output files
SKIF_path = skf.get_SKIF_directory() #get SKIF project root directory
#TablesPath = skf.path_in_project('/' + station + '/TechReports/tabl/')#, your_sys='Mac OC')
#FigPath = skf.path_in_project('/' + station + '/TechReports/pic/')
wfrPath = skf.path_in_project('/' + station + '/fields_' + station + '/')
#Diamond_T_path = skf.path_in_project('/' + station + '/crystals_data_' + station + '/diamond_T/')

wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
#spec1FileName = 'wfr_spec1_' + station + '.wfr' #for spec1
#spec2FileName = 'wfr_spec2_' + station + '.wfr' #for spec2
#stkPFileName = 'stkP.wfr'#for power density

#wfrFileName = [spec1FileName, spec2FileName]#, stkPFileName]
#harmonics number
harm1 = 3

Length = 2.3 # m
undper = 0.018 # m
numper = 128
magf = 1.33

#**********************Output files
SKIF_path = skf.get_SKIF_directory() #get SKIF project root directory
#TablesPath = skf.path_in_project('/' + station + '/TechReports/tabl/')#, your_sys='Mac OC')
#FigPath = skf.path_in_project('/' + station + '/TechReports/pic/')
wfrPath = skf.path_in_project('/' + station + '/fields_' + station + '/')
#Diamond_T_path = skf.path_in_project('/' + station + '/crystals_data_' + station + '/diamond_T/')

wfrPathName = SKIF_path + '/' + station + '/fields_' + station + '/' #example data sub-folder name
wfr1FileName = 'wfr_DCM.wfr' #for DCM

wfrFileName = [wfr1FileName]

#***********Undulator
harmB1 = SRWLMagFldH() #magnetic field harmonic
harmB1.n = 1 #harmonic number
harmB1.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB1.B = magf #magnetic field amplitude [T]

und1 = SRWLMagFldU([harmB1])
und1.per = undper  #period length [m]
und1.nPer = numper #number of periods (will be rounded to integer)

K = 0.9336 * magf * undper * 100 #undulator parameter
E1 = round(4*np.pi*speed_of_light*h_bar*gamma**2/(undper*(1 + K**2/2)), 2) #energy of the first harmonic

print("K = ", K)#, "\n", 'Delta_theta_{} = '.format(11), Delta_theta)

for i in range(1, 25, 2):
    Delta_theta = np.sqrt(4*np.pi*speed_of_light*h_bar/(i*E1)/undper/numper) #angle divergence (the first minimum of intensity)
    print('E{} = '.format(i), round(i*E1, 2), '  ang_{} = '.format(i), round(Delta_theta,7))

magFldCnt = SRWLMagFldC([und1], array('d', [0]), array('d', [0.0000]), array('d', [0.0001])) #Container of all Field Elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = gamma#3./0.51099890221e-03 #relative energy 3 Gev??
sigEperE = 8.6e-04 #relative RMS energy spread
sigX = 33.0e-06 #horizontal RMS size of e-beam [m]
sigXp = 2.65e-06 #horizontal RMS angular divergence [rad]
sigY = 8.6e-07 #vertical RMS size of e-beam [m]
sigYp = 5.0e-07 #vertical RMS angular divergence [rad]
#2nd order stat. moments:
eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2> 
eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2> 
eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

#***********Auxiliary Electron Trajectory structure (for test)
partTraj = SRWLPrtTrj() #defining auxiliary trajectory structure
partTraj.partInitCond = eBeam.partStatMom1
partTraj.allocate(20001) 
partTraj.ctStart = -1.6 #Start "time" for the calculation
partTraj.ctEnd = 1.6

#***********Precision Parameters
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = 21 #final UR harmonic to take into account
arPrecF[2] = 1.5 #longitudinal integration precision parameter
arPrecF[3] = 1.5 #azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.001 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]
#%%

mesh_wfr = 250
distance = 30. #[m]
a = 0.0015 #[m]

wfr1 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr1.allocate(1, mesh_wfr, mesh_wfr) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = distance  #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = round(harm1*E1)#4205 #Initial Photon Energy [eV]
wfr1.mesh.eFin = wfr1.mesh.eStart #Final Photon Energy [eV]
wfr1.mesh.xStart = -a#*distance*1e-6#-a/distance #Initial Horizontal Position [m]
wfr1.mesh.xFin = a#*distance*1e-6#a/distance #Final Horizontal Position [m]
wfr1.mesh.yStart = -a#*distance*1e-6#-a/distance #Initial Vertical Position [m]
wfr1.mesh.yFin = a#*distance*1e-6#a/distance #Final Vertical Position [m]
wfr1.partBeam = eBeam

wfrContainer = [wfr1]

#%%
            #### Electric field calculation #####
for wfr in wfrContainer:
    print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
    srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecPar)
    print('done')

#%%
            ######### Intensity Ploting#######
skf.skf_wfr_subplot_XY(wfr1, fourth_plot=0)


#%% A Long long part of the code with the crystal definition
###   Double Crystal Monochromator (for harm 17)
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
orientDataCr1 = DCM_Cr1.find_orient(wfr1.avgPhotEn)
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
orientDataCr2 = DCM_Cr2.find_orient(wfr1.avgPhotEn, _ang_dif_pl=3.141593)
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

Drift_BEFORE_Cr = SRWLOptD(1.5) #Drift before silicon(111) DCM for 3 harmonic
Drift_AFTER_Cr = SRWLOptD(1.5) #Drift before silicon(111) DCM for 3 harmonic

SSA = SRWLOptA('r', 'a', 2*slit_x*1e-03, 2*slit_y*1e-03) #slit parameter
#%%
#Wavefront Propagation Parameters:
#                       [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
pprPar0 =               [ 0,  0,  1.,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppSSA               =   [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppDrift_BEFORE_Cr  =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppDrift_AFTER_Cr  =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

prParPost =             [ 0,  0,  1.,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

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

optBL = SRWLOptC([DCM_Cr1, DCM_Cr2, Drift_AFTER_Cr],
                 [pprPar0, pprPar0, ppDrift_AFTER_Cr, prParPost])

#%%
#///////////WAVEFRONT PROPOGATION//////#
print('   Simulating Electric Field Wavefront Propagation ... ', end='')
t0 = time.time()
srwl.PropagElecField(wfr1, optBL) # 3 harm for DCM
print('done; lasted', round(time.time() - t0), 's')

skf.skf_wfr_subplot_XY(wfr1, fourth_plot=0)


