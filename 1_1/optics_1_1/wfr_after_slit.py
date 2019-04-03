#############################################################################
#Extract two electric fields from files. Plot intensity distribution and spectrum.
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle

print('SKIF Extended Example for 1-1 # 2:')
print('Extract two electric fields from files. Plot intensity distribution and spectrum.')

#**********************Input Parameters:
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
wfr1FileName = 'spec1_1.wfr' # for spectrum
wfr2FileName = 'wfr_harm2.wfr' #for harm2
wfr3FileName = 'wfr_harm3.wfr' #for harm3
wfr4FileName = 'wfr_harm4.wfr' #for harm4
wfr5FileName = 'wfr_harm5.wfr' #for harm5

#%%
print('extracting wfr from files')
afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr2   =  pickle.load(afile)
afile.close()

#afile = open(wfrPathName + wfr3FileName, 'rb')
#wfr3 =  pickle.load(afile)
#afile.close()
#
#afile = open(wfrPathName + wfr4FileName, 'rb')
#wfr4 =  pickle.load(afile)
#afile.close()
#
#afile = open(wfrPathName + wfr5FileName, 'rb')
#wfr5 =  pickle.load(afile)
#afile.close()
print('finishing extracting')

p_wfr1 = deepcopy(wfr1)
p_wfr2 = deepcopy(wfr2)
#p_wfr3 = deepcopy(wfr3)
#p_wfr4 = deepcopy(wfr4)
#p_wfr5 = deepcopy(wfr5)

print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
uti_plot2d(arI2, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])

arI2x = array('f', [0]*wfr2.mesh.nx) #array to take 1D intensity data (vs X)
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2x, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(horizontal cut at x = 0)'])

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

#%%

#Diamond(111) Crystal Constants:
dSpC400 = 2.0592929401455575#(Diamond(111))#0.89178 #Crystal reflecting planes d-spacing for Diamond(111) crystal
#psi0rC400 = -0.17530e-04; psi0iC400 = 0.21089e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
psi0rC400 = -7.4441e-6; psi0iC400 = 2.9887e-9#-0.21732e-04; psi0iC400 = 0.28005e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
#psihrC400 = -0.45300e-05; psihiC400 = 0.20314E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihrC400 = 2.707e-6;   psihiC400 = 2.0807e-9#-0.54377e-05; psihiC400 = 0.25934E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihbrC400 = psihrC400; psihbiC400 = psihiC400 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
thickCryst = 1000e-06 #0.5e-03 #Thickness of each crystal [m]
angAsCryst = 0 #Asymmetry angle of each crystal [rad]
#1st Crystal:
opCr1 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst, _uc=2)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):

orientDataCr1 = opCr1.find_orient(3050)#wfr1.avgPhotEn+123)
orientCr1 = orientDataCr1[0] #1st crystal orientation
tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
#print('   1st crystal orientation:'); print('   t=', tCr1, 's=', orientCr1[1], 'n=', nCr1)
#Set crystal orientation:
opCr1.set_orient()#nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
print(nCr1[0], nCr1[1], nCr1[2])
orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
uti_math.matr_print(TrM)


#***************** Optical Elements and Propagation Parameters
slit = 0.2 # mm
Drift_BEFORE_SLIT = SRWLOptD(1.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)
SSA = SRWLOptA('r', 'a', 2*slit*1e-03, 2*slit*1e-03)
Drift_AFTER_SLIT = SRWLOptD(1.)
#Drift_SSA_SCREEN = SRWLOptD(1.) #Drift from SSA to Center of Vertically Focusing K-B Mirror (VKB)

#Wavefront Propagation Parameters:
#                       [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
prParInit =             [ 0,  0, 1.,   1,  0, 1.,  2.,  1.3,  2.,  0,  0,   0]
ppDrift_BEFORE_SLIT =   [ 0,  0, 1.0,  1,  0, 1.,  0.7, 1.,   0.7, 0,  0,   0]
ppSSA               =   [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]
ppDrift_AFTER_SLIT  =   [ 0,  0, 1.0,  1,  0, 1.1, 1.0, 1.1,  1.0, 0,  0,   0]
prPar0 =                [ 0,  0,  1.,  1,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]
prParPost =             [ 0,  0,  1.,  0,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]

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
optBL = SRWLOptC([opCr1, Drift_BEFORE_SLIT], 
                 [prParInit, prPar0, ppDrift_BEFORE_SLIT, prParPost])
#optBL = SRWLOptC([Drift_BEFORE_SLIT], [ppDrift_BEFORE_SLIT, prParPost])

#%%
#///////////WAVEFRONT PROPOGATION//////#
print('   Simulating Electric Field Wavefront Propagation ... ', end='')
t0 = time.time();
srwl.PropagElecField(wfr1, optBL)
srwl.PropagElecField(wfr2, optBL)
#srwl.PropagElecField(p_wfr3, optBL)
#srwl.PropagElecField(p_wfr4, optBL)
#srwl.PropagElecField(p_wfr5, optBL)
print('done; lasted', round(time.time() - t0), 's')

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

arI2x = array('f', [0]*wfr2.mesh.nx) #array to take 1D intensity data (vs X)
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2x, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(horizontal cut at x = 0)'])

print('   Extracting Intensity from the Propagated Electric Field  ... ', end='')
arI1 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" 2D array to take intensity data
srwl.CalcIntFromElecField(arI1, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
uti_plot2d(arI1, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])






#%%
'''
            ######### Spectrum #######
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
p_arI1 = array('f', [0]*p_wfr1.mesh.ne)
srwl.CalcIntFromElecField(p_arI1, wfr1, 6, 2, 0, p_wfr1.mesh.eStart, p_wfr1.mesh.xStart, p_wfr1.mesh.yStart)
uti_plot1d(p_arI1, [p_wfr1.mesh.eStart, p_wfr1.mesh.eFin, p_wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

#%%
            ######## Intensity #######
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
uti_plot2d(arI2, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])

p_arI2 = array('f', [0]*p_wfr2.mesh.nx*p_wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(p_arI2, p_wfr2, 6, 0, 3, p_wfr2.mesh.eStart, 0, 0)
uti_plot2d(p_arI2, [1000*p_wfr2.mesh.xStart, 1000*p_wfr2.mesh.xFin, p_wfr2.mesh.nx], [1000*p_wfr2.mesh.yStart, 1000*p_wfr2.mesh.yFin, p_wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(p_wfr2.mesh.eStart) + ' eV'])
print('done')

#print('   Extracting Intensity from calculated Electric Field ... ', end='')
#arI3 = array('f', [0]*wfr3.mesh.nx*wfr3.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(arI3, wfr3, 6, 0, 3, wfr3.mesh.eStart, 0, 0)
#uti_plot2d(arI3, [1000*wfr3.mesh.xStart, 1000*wfr3.mesh.xFin, wfr3.mesh.nx], [1000*wfr3.mesh.yStart, 1000*wfr3.mesh.yFin, wfr3.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr3.mesh.eStart) + ' eV'])
#
#p_arI3 = array('f', [0]*p_wfr3.mesh.nx*p_wfr3.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(p_arI3, p_wfr3, 6, 0, 3, p_wfr3.mesh.eStart, 0, 0)
#uti_plot2d(p_arI3, [1000*p_wfr3.mesh.xStart, 1000*p_wfr3.mesh.xFin, p_wfr3.mesh.nx], [1000*p_wfr3.mesh.yStart, 1000*p_wfr3.mesh.yFin, p_wfr3.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(p_wfr3.mesh.eStart) + ' eV'])
#print('done')
#
#print('   Extracting Intensity from calculated Electric Field ... ', end='')
#arI4 = array('f', [0]*wfr4.mesh.nx*wfr4.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(arI4, wfr4, 6, 0, 3, wfr4.mesh.eStart, 0, 0)
#uti_plot2d(arI4, [1000*wfr4.mesh.xStart, 1000*wfr4.mesh.xFin, wfr4.mesh.nx], [1000*wfr4.mesh.yStart, 1000*wfr4.mesh.yFin, wfr4.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr4.mesh.eStart) + ' eV'])
#
#p_arI4 = array('f', [0]*p_wfr4.mesh.nx*p_wfr4.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(p_arI4, p_wfr4, 6, 0, 3, p_wfr4.mesh.eStart, 0, 0)
#uti_plot2d(p_arI4, [1000*p_wfr4.mesh.xStart, 1000*p_wfr4.mesh.xFin, p_wfr4.mesh.nx], [1000*p_wfr4.mesh.yStart, 1000*p_wfr4.mesh.yFin, p_wfr4.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(p_wfr4.mesh.eStart) + ' eV'])
#print('done')
#
#print('   Extracting Intensity from calculated Electric Field ... ', end='')
#arI5 = array('f', [0]*wfr5.mesh.nx*wfr5.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(arI5, wfr5, 6, 0, 3, wfr5.mesh.eStart, 0, 0)
#uti_plot2d(arI5, [1000*wfr5.mesh.xStart, 1000*wfr5.mesh.xFin, wfr5.mesh.nx], [1000*wfr5.mesh.yStart, 1000*wfr5.mesh.yFin, wfr5.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr5.mesh.eStart) + ' eV'])
#
#p_arI5 = array('f', [0]*p_wfr5.mesh.nx*p_wfr5.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(p_arI5, p_wfr5, 6, 0, 3, p_wfr5.mesh.eStart, 0, 0)
#uti_plot2d(p_arI5, [1000*p_wfr5.mesh.xStart, 1000*p_wfr5.mesh.xFin, p_wfr5.mesh.nx], [1000*p_wfr5.mesh.yStart, 1000*p_wfr5.mesh.yFin, p_wfr5.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(p_wfr5.mesh.eStart) + ' eV'])
#print('done')



'''







