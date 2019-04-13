#############################################################################
#Create a tapered undulator/wiggler magnetic structure. Calculate and draw !two! electric field files 
#otimised for extracting spectrum and intensity and the trajectory of the electron. Save it to files using pickle lib
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
print('SRWLIB Extended Example # 6:')
print('Create a tapered undulator/wiggler magnetic structure. Calculate and draw !two! electric field files otimised for extracting spectrum and intensity and the trajectory of the electron. Save it to files using pickle lib')

#**********************Output files
PathName = '/home/andrei/Documents/SKIF_XAS_beamline/fields/' #example data sub-folder name
FileName = 'undulator_wiggler_traj.trj' #file name for output UR flux data

#***********Undulator
undarr = []
undarrH = []
distz =  []
distx =  []
disty =  []

Length = 2.3 # m
undper = 0.0186 # m
numper = 30
magf = 1
magf_step = 1/100

PER = 140
NumPIECE = 10
undper = 0.017
numper = PER/NumPIECE
magf = 1
magf_step = 1.8/100

Bx = 1 #Peak Horizontal field [T]
By = 1.15 #+ magf_step #+ (NumPIECE-1)*magf_step#Peak Vertical field [T]
phBx = 0#Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx = 1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy = 1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = -0.00 #Transverse Coordinates of Undulator Center [m]
ycID = -0.0
zcID = 0 #Longitudinal Coordinate of Undulator Center [m]

K = 0.965 * magf * undper * 100
print("K = ",K)

for i in range(NumPIECE):
    B2 = By + i*magf_step*By
    #phBx = rn.uniform(0, 2*np.pi)
    und = SRWLMagFldU([SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undper, numper)#, SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
    undarr.append(und)
    distz.append((i)*(undper*numper + 3*undper))
    distx.append(0)
    disty.append(0)
print(distx, distz)

magFldCnt = SRWLMagFldC(undarr, array('d', distx), array('d', disty), array('d', distz)) #Container of all Field Elements

'''
harmB1 = SRWLMagFldH() #magnetic field harmonic
#harmB1.n = 0 #harmonic number
harmB1.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB1.B = magf #magnetic field amplitude [T]
und1 = SRWLMagFldU([harmB1])
und1.per = undper  #period length [m]
und1.nPer = numper #number of periods (will be rounded to integer)

harmB2 = SRWLMagFldH() #magnetic field harmonic
harmB2.n = 1 #harmonic number
harmB2.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB2.B = magf + 1*magf_step #magnetic field amplitude [T]
harmB2.ph = rn.uniform(0, 2*np.pi)
und2 = SRWLMagFldU([harmB2])
und2.per = undper  #period length [m]
und2.nPer = numper #number of periods (will be rounded to integer)

harmB3 = SRWLMagFldH() #magnetic field harmonic
harmB3.n = 1 #harmonic number
harmB3.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB3.B = magf + 2*magf_step #magnetic field amplitude [T]
harmB3.ph = rn.uniform(0, 2*np.pi)
und3 = SRWLMagFldU([harmB3])
und3.per = undper  #period length [m]
und3.nPer = numper #number of periods (will be rounded to integer)

harmB4 = SRWLMagFldH() #magnetic field harmonic
harmB4.n = 1 #harmonic number
harmB4.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB4.B = magf + 3*magf_step #magnetic field amplitude [T]
harmB4.ph = rn.uniform(0, 2*np.pi)
und4 = SRWLMagFldU([harmB4])
und4.per = undper  #period length [m]
und4.nPer = numper #number of periods (will be rounded to integer)

harmB5 = SRWLMagFldH() #magnetic field harmonic
harmB5.n = 1 #harmonic number
harmB5.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB5.B = magf + 4*magf_step #magnetic field amplitude [T]
harmB5.ph = rn.uniform(0, 2*np.pi)
und5 = SRWLMagFldU([harmB5])
und5.per = undper  #period length [m]
und5.nPer = numper #number of periods (will be rounded to integer)
'''

#magFldCnt = SRWLMagFldC([und1, und2, und3, und4, und5], 
#                        array('d', [0, 0, 0, 0, 0]), 
#                        array('d', [0, 0, 0, 0, 0]), 
#                        array('d', [0, numper*undper+2*undper*rn.uniform(1, 2), 2*numper*undper+4*undper*rn.uniform(1, 2), 3*numper*undper+6*undper*rn.uniform(1, 2), 4*numper*undper+8*undper*rn.uniform(1, 2)])) #Container of all Field Elements
#magFldCnt = SRWLMagFldC([und1, und2, und3, und4, und5], 
#                        array('d', [0, 0, 0, 0, 0]), 
#                        array('d', [0, 0, 0, 0, 0]), 
#                        array('d', [0, numper*undper+3*undper, 2*numper*undper+6*undper, 3*numper*undper+9*undper, 4*numper*undper+12*undper])) #Container of all Field Elements

#magFldCnt = SRWLMagFldC([und1, und2, und3], array('d', [0, 0, 0]), array('d', [0, 0, 0]), array('d', [0, numper*undper+2*undper*rn.uniform(1, 2), 2*numper*undper+4*undper*rn.uniform(1, 2)])) #Container of all Field Elements
#magFldCnt = SRWLMagFldC([und1, und2], array('d', [0, 0]), array('d', [0, 0]), array('d', [0, numper*undper+3*undper])) #Container of all Field Elements
#magFldCnt = SRWLMagFldC([und1], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all Field Elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = 3./0.51099890221e-03 #relative energy 3 Gev??
sigEperE = 0.001#0.00089 #relative RMS energy spread
sigX = 33.33e-06 #horizontal RMS size of e-beam [m]
sigXp = 16.5e-06 #horizontal RMS angular divergence [rad]
sigY = 2.912e-06 #vertical RMS size of e-beam [m]
sigYp = 2.7472e-06 #vertical RMS angular divergence [rad]
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

#***********UR Stokes Parameters (mesh) for Spectral Flux

wfr1 = SRWLWfr() #For spectrum vs photon energy
wfr1.allocate(550, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 4000 #Initial Photon Energy [eV]
wfr1.mesh.eFin = 10500#4300. #Final Photon Energy [eV]
#wfr1.avgPhotEn= #4205
a = 0.0009
wfr1.mesh.xStart = -a #Initial Horizontal Position [m]
wfr1.mesh.xFin = a #Final Horizontal Position [m]
wfr1.mesh.yStart = -a #Initial Vertical Position [m]
wfr1.mesh.yFin = a #Final Vertical Position [m]
wfr1.partBeam = eBeam

wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr2.allocate(1, 151, 151) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = 1400#4205 #Initial Photon Energy [eV]
wfr2.mesh.eFin = wfr2.mesh.eFin #Final Photon Energy [eV]
a = 0.001
wfr2.mesh.xStart = -a #Initial Horizontal Position [m]
wfr2.mesh.xFin = a #Final Horizontal Position [m]
wfr2.mesh.yStart = -a #Initial Vertical Position [m]
wfr2.mesh.yFin = a #Final Vertical Position [m]
wfr2.partBeam = eBeam
print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
print('done')


#print('   Performing Electric Field (wavefront at fixed photon energy) calculation ... ', end='')
#srwl.CalcElecFieldSR(wfr2, 0, magFldCnt, arPrecPar)
#print('done')

            ######### Spectrum #######
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
print('done')


print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')

#uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
E = np.linspace(wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne)
plt.figure(figsize=(1.5*4,1.5*3))
plt.plot(E, arI1, color='blue')
plt.xlabel(r'$E, [эВ]$', fontsize=14, labelpad = 0.0)
y = plt.ylabel(r'$I, [\gamma/с/0.1\%пп/мм^{2}]$', fontsize=14, labelpad = 0.0, rotation=90)
plt.title('')
#y.set_rotation(0)
ax = plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['top'].set_position(('axes',0))
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_position(('data',0))
ax.xaxis.set_label_coords(0.95, -0.08)
ax.yaxis.set_label_coords(-0.05, 0.7)
#plt.xticks([1200, 1300, 1400, 1500, 1600, 1700, 1800],fontsize=12)
        #  [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$'])
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
#plt.yticks([1e14,2e14,3e14,4e14,5e14,6e14,7e14], fontsize=12)
#plt.ylim(0, 7e14)
plt.xlim(wfr1.mesh.eStart, wfr1.mesh.eFin)
plt.savefig('/home/andrei/Documents/SKIF_XAS_beamline/TeXDoc/pic/' +'spec_SRW'+'.pdf')
plt.show()

#            ######### Intensity #######
#print('   Extracting Intensity from calculated Electric Field ... ', end='')
#arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
#srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
#print('done')
#uti_plot2d(arI2, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])

#%%
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
partTraj.ctStart = -0.55 #Start Time for the calculation
partTraj.ctEnd = 3.15#magFldCnt.arMagFld[0].rz

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

#%%
#*****************Saving to files
afile = open(PathName + FileName, 'wb')
pickle.dump(partTraj, afile)
afile.close()
















