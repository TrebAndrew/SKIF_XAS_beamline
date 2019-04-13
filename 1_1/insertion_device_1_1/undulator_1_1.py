#############################################################################
#Stastation 1-1 $$ insertion device
#Create a undulator magnetic structure. Calculate and save wave fronts at desired harmonics. 
#Calculate spectrum
#v0.1

#harmonic numbers 11-th, 13-th, 17-th and 23-th
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

print('SKIF Extended Example for 1-1 # 1:')
print('Create an undulator for 1-1 station.')
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
#**********************Output files
PathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
FileName = 'undulator_traj.trj' #file name for output electrom traj data

wfrPathName = '/home/andrei/Documents/SKIF_XAS_beamline/1_1/fields_1_1/' #example data sub-folder name
wfr1FileName = 'spec1_1_i.wfr' #for spectrum
wfr2FileName = 'wfr_harm2.wfr' #for harm2
wfr3FileName = 'wfr_harm3.wfr' #for harm3
wfr4FileName = 'wfr_harm4.wfr' #for harm4
wfr5FileName = 'wfr_harm5.wfr' #for harm5
wfr6FileName = 'spec1_1_ii.wfr' #for harm5
stkPFileName = 'stkP_harm5.wfr'#for power dens

wfrFileName = [wfr1FileName, wfr2FileName, wfr3FileName, wfr4FileName, wfr5FileName, wfr6FileName]#, stkPFileName]

#***********Undulator
harmB1 = SRWLMagFldH() #magnetic field harmonic
harmB1.n = 1 #harmonic number
harmB1.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB1.B = magf #magnetic field amplitude [T]

und1 = SRWLMagFldU([harmB1])
und1.per = undper  #period length [m]
und1.nPer = numper #number of periods (will be rounded to integer)

K = 0.9336 * magf * undper * 100
E1 = round(4*np.pi*speed_of_light*h_bar*gamma**2/(undper*(1 + K**2/2)), 2)

print("K = ", K)#, "\n", 'Delta_theta_{} = '.format(11), Delta_theta)

for i in range(1, 25, 2):
    Delta_theta = np.sqrt(4*np.pi*speed_of_light*h_bar/(i*E1)/undper/numper)
    print('E{} = '.format(i), round(i*E1, 2), '  ang_{} = '.format(i), round(Delta_theta,7))

magFldCnt = SRWLMagFldC([und1], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all Field Elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = gamma#3./0.51099890221e-03 #relative energy 3 Gev??
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
#%%
mesh_wfr = 150
distance = 20.
a = 14
#***********UR Stokes Parameters (mesh) for Spectral Flux
wfr1 = SRWLWfr() #For spectrum vs photon energy
wfr1.allocate(60000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = distance #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 100 #Initial Photon Energy [eV]
wfr1.mesh.eFin = 20000#4300. #Final Photon Energy [eV]
wfr1.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr1.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr1.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr1.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr1.partBeam = eBeam

afile = open(wfrPathName + 'spec.wfr', 'wb')
pickle.dump(wfr1, afile)
afile.close()

wfr6 = SRWLWfr() #For spectrum vs photon energy
wfr6.allocate(2000, 25, 25) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr6.mesh.zStart = distance #Longitudinal Position [m] at which SR has to be calculated
wfr6.mesh.eStart = 7000 #Initial Photon Energy [eV]
wfr6.mesh.eFin = 15000#4300. #Final Photon Energy [eV]
wfr6.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr6.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr6.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr6.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr6.partBeam = eBeam

#a = 14#*distance*1e-6 #m
wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr2.allocate(1, mesh_wfr, mesh_wfr) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = distance  #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = round(harm2*E1)#4205 #Initial Photon Energy [eV]
wfr2.mesh.eFin = wfr2.mesh.eStart #Final Photon Energy [eV]
wfr2.mesh.xStart = -a*distance*1e-6#-a/distance #Initial Horizontal Position [m]
wfr2.mesh.xFin = a*distance*1e-6#a/distance #Final Horizontal Position [m]
wfr2.mesh.yStart = -a*distance*1e-6#-a/distance #Initial Vertical Position [m]
wfr2.mesh.yFin = a*distance*1e-6#a/distance #Final Vertical Position [m]
wfr2.partBeam = eBeam

wfr3 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr3.allocate(1, mesh_wfr, mesh_wfr) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr3.mesh.zStart = distance  #Longitudinal Position [m] at which SR has to be calculated
wfr3.mesh.eStart = round(harm3*E1)#4205 #Initial Photon Energy [eV]
wfr3.mesh.eFin = wfr3.mesh.eStart #Final Photon Energy [eV]
wfr3.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr3.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr3.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr3.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr3.partBeam = eBeam

wfr4 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr4.allocate(1, mesh_wfr, mesh_wfr) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr4.mesh.zStart = distance  #Longitudinal Position [m] at which SR has to be calculated
wfr4.mesh.eStart = round(harm4*E1)#4205 #Initial Photon Energy [eV]
wfr4.mesh.eFin = wfr4.mesh.eStart #Final Photon Energy [eV]
wfr4.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr4.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr4.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr4.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr4.partBeam = eBeam

wfr5 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr5.allocate(1, mesh_wfr, mesh_wfr) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr5.mesh.zStart = distance  #Longitudinal Position [m] at which SR has to be calculated
wfr5.mesh.eStart = round(harm5*E1)#4205 #Initial Photon Energy [eV]
wfr5.mesh.eFin = wfr5.mesh.eStart #Final Photon Energy [eV]
wfr5.mesh.xStart = -a*distance*1e-6 #Initial Horizontal Position [m]
wfr5.mesh.xFin = a*distance*1e-6 #Final Horizontal Position [m]
wfr5.mesh.yStart = -a*distance*1e-6 #Initial Vertical Position [m]
wfr5.mesh.yFin = a*distance*1e-6 #Final Vertical Position [m]
wfr5.partBeam = eBeam

stkP = SRWLStokes() #for power density
stkP.allocate(1, 101, 101) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = distance #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -0.02 #initial horizontal position [m]
stkP.mesh.xFin = 0.02 #final horizontal position [m]
stkP.mesh.yStart = -0.015 #initial vertical position [m]
stkP.mesh.yFin = 0.015 #final vertical position [m]

wfrContainer = [wfr1, wfr2, wfr3, wfr4, wfr5, wfr6]#, stkP]

#%%
somelist = wfrContainer
somelist = [x for x in somelist if x not in [stkP, wfr6]]
#%%
            #### Electric field calculation #####
for wfr in somelist:
    print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
    srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecPar)
    print('done')


#%%
afile = open(wfrPathName + 'spec.wfr', 'wb')
pickle.dump(wfr1, afile)
afile.close()

            ######### Spectrum Ploting#######
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
print('done')

print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
#%%
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI6 = array('f', [0]*wfr6.mesh.ne)
srwl.CalcIntFromElecField(arI6, wfr6, 6, 2, 0, wfr6.mesh.eStart, wfr6.mesh.xStart, wfr6.mesh.yStart)
print('done')

print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
uti_plot1d(arI6, [wfr6.mesh.eStart, wfr6.mesh.eFin, wfr6.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

#%%
            ######### Intensity Ploting#######
print('   Extracting Intensity from calculated Electric Field ... ', end='')
A = 1e6/distance

arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)

uti_plot2d(arI2, [A*wfr2.mesh.xStart, A*wfr2.mesh.xFin, wfr2.mesh.nx], [A*wfr2.mesh.yStart, A*wfr2.mesh.yFin, wfr2.mesh.ny], [r'$Horizontal Position [\mu rad]$', r'$Vertical Position [\mu rad]$', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])
print('done')

arI2x = array('f', [0]*wfr2.mesh.nx) #array to take 1D intensity data (vs X)
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, wfr2.mesh.eStart, 0, 0)
x = np.linspace(A*wfr2.mesh.xStart, A*wfr2.mesh.xFin, wfr2.mesh.nx)

#%%
            ######## Power ########
a = 14
stkP = SRWLStokes() #for power density
stkP.allocate(1, 201, 201) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = distance #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -a*distance*1e-6#-0.02 #initial horizontal position [m]
stkP.mesh.xFin = a*distance*1e-6#0.02 #final horizontal position [m]
stkP.mesh.yStart = -a*distance*1e-6#-0.015 #initial vertical position [m]
stkP.mesh.yFin = a*distance*1e-6#0.015 #final vertical position [m]


print('   Performing Power Density calculation (from field) ... ', end='')
srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
#srwl.CalcPowDenSR(stkPx, eBeam, 0, magFldCnt, arPrecP)
#srwl.CalcPowDenSR(stkPy, eBeam, 0, magFldCnt, arPrecP)
print('done')

afile = open(wfrPathName + stkPFileName, 'wb')
pickle.dump(stkP, afile)
afile.close()

#%% 
print('saving to the files')
#*****************Saving to files
for (wfr, fname) in zip(wfrContainer, wfrFileName):
    print(wfr, fname, "\n")
    afile = open(wfrPathName + fname, 'wb')
    pickle.dump(wfr, afile)
    afile.close()

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
partTraj.ctStart = -1.45 #Start Time for the calculation
partTraj.ctEnd = 1.35#magFldCnt.arMagFld[0].rz

#**********************Calculation traj (SRWLIB function call)
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

'''
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
'''












