#############################################################################
#Create a tapered undulator magnetic structure. Calculate !two! electric field files 
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
print('SKIF Extended Example # 2:')
print('Creating a tapered undulator magnetic structure. Calculate !two! electric field files otimised for extracting spectrum and intensity. Save it to files using pickle lib')

#**********************Output files
wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_4/fields_1_4/' #example data sub-folder name
spec1FileName = 'spec_sectional_1_4.wfr' #file name for output UR flux data
stkPFileName = 'stkP_sectional_1_4.wfr'#for power densitys

wfrFileName = [spec1FileName]#, stkPFileName]

#***********Undulator
undarr = []
distz =  []
distx =  []
disty =  []

PER = 150
NumPIECE = 5
undper = 0.02 # [m]
numper = PER/NumPIECE
magf = 1.0 #[T]
magf_step = (1.5/100)*magf

Bx = 1 #Peak Horizontal field [T]
By = magf # + magf_step #+ (NumPIECE-1)*magf_step#Peak Vertical field [T]
phBx = 0#Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx = 1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy = 1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = -0.00 #Transverse Coordinates of Undulator Center [m]
ycID = -0.0
zcID = 0 #Longitudinal Coordinate of Undulator Center [m]

for i in range(NumPIECE-7, NumPIECE-2):
    B2 = By + (i-1)*magf_step
    print(B2)
    #phBx = rn.uniform(0, 2*np.pi)
    und = SRWLMagFldU([SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undper, numper)#, SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
    undarr.append(und)
    distz.append((i)*(undper*numper + 4*1*undper))#*numper*rn.uniform(0.5, 1.5)) )
    distx.append(0)
    disty.append(0)
print(distx, distz)

K = 0.965 * magf * undper * 100
print("K = ", round(K))
print("Undulator Length = ", PER * undper)
magFldCnt = SRWLMagFldC(undarr, array('d', distx), array('d', disty), array('d', distz)) #Container of all Field Elements
#magFldCnt = SRWLMagFldC([und1, und2, und3], array('d', [0, 0, 0]), array('d', [0, 0, 0]), array('d', [0, numper*undper+2*undper, 2*numper*undper+4*undper])) #Container of all Field Elements
#magFldCnt = SRWLMagFldC([und1], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all Field Elements
#***********

#%%
#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4          #average current [A]
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


wfr1 = SRWLWfr() #For spectrum vs photon energy (Spectral Flux)
wfr1.allocate(4000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = 25. #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 4300. #Initial Photon Energy [eV]
wfr1.mesh.eFin = 5200. #Final Photon Energy [eV]
a = 0.0002
wfr1.mesh.xStart = -a #Initial Horizontal Position [m]
wfr1.mesh.xFin = a #Final Horizontal Position [m]
wfr1.mesh.yStart = -a #Initial Vertical Position [m]
wfr1.mesh.yFin = a #Final Vertical Position [m]
wfr1.partBeam = eBeam

#%%
#**********************Calculation (SRWLIB function calls)
            ######### Spectrum #######
print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
print('done')

E, spec = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV', emittance=0)
skf.skf_plot(E, spec, elec_fld_units='W/mm^2/eV', color='blue', grid=True, 
             linewidth=1, save_fig=True, figure_name='sim_und_spec.pdf', file_path='/home/andrei/Documents/diploma/TexPresent/pic/')
plt.show()
##%% 
#*****************Saving to files
afile = open(wfrPathName + spec1FileName, 'wb')
pickle.dump(wfr1, afile)
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
'''
print('   Performing Power Density calculation (from field) ... ', end='')
srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
print('done')

skf.skf_power_subplot_XY(stkP, units='urad')

afile = open(wfrPathName + stkPFileName, 'wb')
pickle.dump(stkP, afile)
afile.close()
'''

#%% ###### Trajectory #######

numPer = 40 #Number of ID Periods (without counting for terminations)
xcID = 0 #Transverse Coordinates of ID Center [m]
ycID = 0
zcID = 0 #Longitudinal Coordinate of ID Center [m]

part = SRWLParticle()
part.x = 0.00 #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
part.y = 0.000
part.xp = 0 #Initial Transverse Velocities
part.yp = 0
part.gamma = 3/0.51099890221e-03 #Relative Energy
part.relE0 = 1 #Electron Rest Mass
part.nq = -1 #Electron Charge

npTraj = 1001 #Number of Points for Trajectory calculation
fieldInterpMeth = 4 #2 #Magnetic Field Interpolation Method, to be entered into 3D field structures below (to be used e.g. for trajectory calculation):
#1- bi-linear (3D), 2- bi-quadratic (3D), 3- bi-cubic (3D), 4- 1D cubic spline (longitudinal) + 2D bi-cubic
arPrecPar = [1] #General Precision parameters for Trajectory calculation:

#**********************Trajectory structure, where the results will be stored
partTraj = SRWLPrtTrj()
partTraj.partInitCond = part
partTraj.allocate(npTraj, True)
partTraj.ctStart = -1.9 #Start Time for the calculation
partTraj.ctEnd = 1.9#magFldCnt.arMagFld[0].rz

#**********************Calculation (SRWLIB function call)
print('   Performing calculation ... ', end='')
partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, arPrecPar)
print('done')

#**********************Plotting results
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
ctMesh = [partTraj.ctStart, partTraj.ctEnd, partTraj.np]
for i in range(partTraj.np):
    partTraj.arX[i] *= 1000
    partTraj.arY[i] *= 1000
    
uti_plot1d(partTraj.arX, ctMesh, ['ct [m]', 'Horizontal Position [mm]'])
uti_plot1d(partTraj.arY, ctMesh, ['ct [m]', 'Vertical Position [mm]'])
