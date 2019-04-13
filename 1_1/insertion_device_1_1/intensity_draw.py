#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 11:03:43 2019

@author: andrei
"""

#############################################################################
#Extract electric fields from files. Plot intensity distribution and spectrum.
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle
import numpy as np
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
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/spec_1_1/' #example data sub-folder name
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/intens_1_1/' #example data sub-folder name

#%%
for filename in os.listdir(specPathName):
    afile = open(specPathName + filename, 'rb')
    wfr1 =  pickle.load(afile)
    print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
    arI1 = array('f', [0]*wfr1.mesh.ne)
    srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
    uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
    print('done')
    afile.close()



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
