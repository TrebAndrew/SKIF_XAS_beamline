#############################################################################
#Create an undulator structure. Calculate !two! electric field files 
#otimised for extracting spectrum and intensity. Save it to files using pickle lib
#v0.1
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

print('SKIF Extended Example # 1:')
print('Create an undulator structure. Calculate !two! electric field files otimised for extracting spectrum and intensity. Save it to files using pickle lib')
#**********************Output files
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_4/fields_1_4/' #example data sub-folder name
spec1FileName = 'wfr_spec1_1_1.wfr' #for spec1
spec2FileName = 'wfr_spec2_1_1.wfr' #for spec2
stkPFileName = 'stkP.wfr'#for power density

wfrFileName = [spec1FileName, spec2FileName]#, stkPFileName]


#***********Undulator
undarr = []
undarrH = []
distz =  []
distx =  []
disty =  []

Length = 2.3 # m
undper = 0.018 # m
numper = 128
magf = 1.3

harmB1 = SRWLMagFldH() #magnetic field harmonic
harmB1.n = 1 #harmonic number
harmB1.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB1.B = magf #magnetic field amplitude [T]

und1 = SRWLMagFldU([harmB1])
und1.per = undper  #period length [m]
und1.nPer = numper #number of periods (will be rounded to integer)

K = 0.965 * magf * undper * 100
print("K = ", round(K))
print("Undulator Length = ", undper * numper)
magFldCnt = SRWLMagFldC([und1], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all Field Elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = 3./0.51099890221e-03 #relative energy 3 Gev??
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
mesh_wfr = 200
distance = 25.
a = 0.001

#***********UR Stokes Parameters (mesh) for Spectral Flux
wfr1 = SRWLWfr() #For spectrum vs photon energy
wfr1.allocate(60000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = distance #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 100 #Initial Photon Energy [eV]
wfr1.mesh.eFin = 60000#4300. #Final Photon Energy [eV]
wfr1.mesh.xStart = -a#*distance*1e-6 #Initial Horizontal Position [m]
wfr1.mesh.xFin = a#*distance*1e-6 #Final Horizontal Position [m]
wfr1.mesh.yStart = -a#*distance*1e-6 #Initial Vertical Position [m]
wfr1.mesh.yFin = a#*distance*1e-6 #Final Vertical Position [m]
wfr1.partBeam = eBeam

#***********UR Stokes Parameters (mesh) for Spectral Flux
wfr2 = SRWLWfr() #For spectrum vs photon energy
wfr2.allocate(20000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = distance #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = 100 #Initial Photon Energy [eV]
wfr2.mesh.eFin = 20000#4300. #Final Photon Energy [eV]
wfr2.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr2.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr2.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr2.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr2.partBeam = eBeam

wfrContainer = [wfr1, wfr2]#, stkP]

#%%exclude unnecessary objects from wfrContainer
#somelist = wfrContainer
#somelist = [x for x in somelist if x not in [stkP]]
#print(somelist)
#%%
            #### Electric field calculation #####
for wfr in wfrContainer:
    print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
    srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecPar)
    print('done')
#%%
            ######### Spectrum Ploting#######
print('   Performing extraction of the wfrs and plotting it ... ', end='')

z1 = wfr1.mesh.zStart
z2 = wfr2.mesh.zStart 

plt.figure()
E, spec1 = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV', emittance=0)
skf.skf_plot(E, spec1, elec_fld_units='W/mm^2/eV', color='blue', grid=True, linewidth=1.5, show=True)
print('Total power density = ', round((np.sum(spec1))*((wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne)), '[W/mm^2]')

plt.figure()
E, spec2 = skf.renorm_wfr(wfr2, elec_fld_units='W/mm^2/eV', emittance=0)
skf.skf_plot(E, spec2, elec_fld_units='W/mm^2/eV', color='blue', grid=True, linewidth=1, show=True)

print('done')
#%% 
print('saving to the files')
#*****************Saving to files
for (wfr, fname) in zip(wfrContainer, wfrFileName):
    print(wfr, fname, "\n")
    afile = open(wfrPathName + fname, 'wb')
    pickle.dump(wfr, afile)
    afile.close()

#%%
######## Power ########

a = 10#[mm]
b = 10#[mm]
A = 1e-3#[mm] -> [m]
distance = 25#[m]

stkP = SRWLStokes() #for power density
stkP.allocate(1, 151, 151) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = distance #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -a*A#*distance*1e-6#-0.02 #initial horizontal position [m]
stkP.mesh.xFin = a*A#*distance*1e-6#0.02 #final horizontal position [m]
stkP.mesh.yStart = -b*A#*distance*1e-6#-0.015 #initial vertical position [m]
stkP.mesh.yFin = b*A#*distance*1e-6#0.015 #final vertical position [m]

print('   Performing Power Density calculation (from field) ... ', end='')
srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
print('done')

skf.skf_power_subplot_XY(stkP, units='urad')

afile = open(wfrPathName + stkPFileName, 'wb')
pickle.dump(stkP, afile)
afile.close()









