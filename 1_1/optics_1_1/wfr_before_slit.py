#############################################################################
#Extract two electric fields from files. Plot intensity distribution and spectrum.
#############################################################################

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
print('Extract two electric fields from files. Plot intensity distribution and spectrum.')

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3.0/0.51099890221e-03
e_ = 1.60218e-19 #elementary charge

distance = 20.
a = 10
#**********************Input Parameters:
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
wfr1FileName = 'spec1_1_i.wfr'   #for spectrum
wfr2FileName = 'wfr_harm2.wfr' #for harm2
wfr3FileName = 'wfr_harm3.wfr' #for harm3
wfr4FileName = 'wfr_harm4.wfr' #for harm4
wfr5FileName = 'wfr_harm5.wfr' #for harm5
wfr6FileName = 'spec1_1_ii.wfr' #for harm5
stkPFileName = 'stkP_harm5.wfr'#for power dens

#%%
print('extracting wfr from files')
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

afile = open(wfrPathName + wfr5FileName, 'rb')
wfr5 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr6FileName, 'rb')
wfr6 =  pickle.load(afile)
afile.close()
#
afile = open(wfrPathName + stkPFileName, 'rb')
stkP =  pickle.load(afile)
afile.close()
print('finishing extracting')

            ######### Spectrum #######
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
print('done')

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI6 = array('f', [0]*wfr6.mesh.ne)
srwl.CalcIntFromElecField(arI6, wfr6, 6, 2, 0, wfr6.mesh.eStart, wfr6.mesh.xStart, wfr6.mesh.yStart)
uti_plot1d(arI1, [wfr6.mesh.eStart, wfr6.mesh.eFin, wfr6.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
print('done')
#%%
            ######### Intensity #######
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
A = 1e6/distance
uti_plot2d(arI2, [A*wfr2.mesh.xStart, A*wfr2.mesh.xFin, wfr2.mesh.nx], [A*wfr2.mesh.yStart, A*wfr2.mesh.yFin, wfr2.mesh.ny], [r'$Horizontal Position [\mu rad]$', r'$Vertical Position [\mu rad]$', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])
    

arI2x = array('f', [0]*wfr2.mesh.nx) #array to take 1D intensity data (vs X)
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2x, [A*wfr2.mesh.xStart, A*wfr2.mesh.xFin, wfr2.mesh.nx], ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(horizontal cut at x = 0)'])

arI2y = array('f', [0]*wfr2.mesh.ny) #array to take 1D intensity data (vs Y)
srwl.CalcIntFromElecField(arI2y, wfr2, 6, 0, 2, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2y, [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Vertical Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(vertical cut at y = 0)'])
    
x = np.linspace(A*wfr2.mesh.xStart, A*wfr2.mesh.xFin, wfr2.mesh.nx)
sigma = skf.calc_bandwidth(x, arI2x)

uti_plot_show() #show all graphs (and block execution)
print('done')

'''
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI3 = array('f', [0]*wfr3.mesh.nx*wfr3.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI3, wfr3, 6, 0, 3, wfr3.mesh.eStart, 0, 0)
uti_plot2d(arI3, [1000*wfr3.mesh.xStart, 1000*wfr3.mesh.xFin, wfr3.mesh.nx], [1000*wfr3.mesh.yStart, 1000*wfr3.mesh.yFin, wfr3.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr3.mesh.eStart) + ' eV'])
print('done')

print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI4 = array('f', [0]*wfr4.mesh.nx*wfr4.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI4, wfr4, 6, 0, 3, wfr4.mesh.eStart, 0, 0)
uti_plot2d(arI4, [1000*wfr4.mesh.xStart, 1000*wfr4.mesh.xFin, wfr4.mesh.nx], [1000*wfr4.mesh.yStart, 1000*wfr4.mesh.yFin, wfr4.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr4.mesh.eStart) + ' eV'])
print('done')

print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI5 = array('f', [0]*wfr5.mesh.nx*wfr5.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI5, wfr5, 6, 0, 3, wfr5.mesh.eStart, 0, 0)
uti_plot2d(arI5, [1000*wfr5.mesh.xStart, 1000*wfr5.mesh.xFin, wfr5.mesh.nx], [1000*wfr5.mesh.yStart, 1000*wfr5.mesh.yFin, wfr5.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr5.mesh.eStart) + ' eV'])
print('done')
'''
#%% 
    
       ######### Power Intensity ###########
            
A = 1e3#1e6/distance
plotMeshX = [A*stkP.mesh.xStart, A*stkP.mesh.xFin, stkP.mesh.nx]
plotMeshY = [A*stkP.mesh.yStart, A*stkP.mesh.yFin, stkP.mesh.ny]

uti_plot2d(stkP.arS, plotMeshX, plotMeshY, [r'$Horizontal Position [\mu rad]$', r'$Vertical Position [\mu rad]$', 'Power Density'])
print(np.sum(stkP.arS)*(1e3*(stkP.mesh.xFin - stkP.mesh.xStart )/stkP.mesh.nx)**2)            

powDenVsX = array('f', [0]*stkP.mesh.nx)
for i in range(stkP.mesh.nx): powDenVsX[i] = stkP.arS[stkP.mesh.nx*int(stkP.mesh.ny*0.5) + i]
uti_plot1d(powDenVsX, plotMeshX, [r'$Horizontal Position [\mu rad]$', 'Power Density [W/mm^2]', 'Power Density\n(horizontal cut at y = 0)'])

powDenVsY = array('f', [0]*stkP.mesh.ny)
for i in range(stkP.mesh.ny): powDenVsY[i] = stkP.arS[int(stkP.mesh.nx*0.5) + i*stkP.mesh.ny]
uti_plot1d(powDenVsY, plotMeshY, [r'$Vertical Position [\mu rad]$', 'Power Density [W/mm^2]', 'Power Density\n(vertical cut at x = 0)'])

uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')
#%%
            ######### Power Intensity Spectrum ###########
afile = open(wfrPathName + 'spec.wfr', 'rb')
wfr1 =  pickle.load(afile)
afile.close()

#%%
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
E1 = np.linspace(wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne)

E, spec = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV')

plt.figure(figsize=(1.5*10,1.5*3))
skf.skf_plot(E, spec,  elec_fld_units='W/mm^2/eV')

T = skf.pycry_trans(crystal='diamond', Emin=wfr1.mesh.eStart, Emax=wfr1.mesh.eFin, ne=wfr1.mesh.ne)
spec = spec[:len(T)]
E = E[:len(T)]
skf.skf_plot(E, spec*T, color='red', elec_fld_units='W/mm^2/eV', grid=True)

print('absorbed power = ', round(np.sum(spec) - np.sum(spec*T)), 'W')









