#############################################################################
#Create a tapered undulator/any magnetic structure. Calculate !two! electric field files 
#otimised for extracting spectrum and intensity. Save it to files using pickle lib
#v0.1
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import pickle

print('SRWLIB Extended Example # 6:')
print('Calculating spectral flux of undulator radiation by finite-emittance electron beam collected through a finite aperture and power density distribution of this radiation (integrated over all photon energies)')

#**********************Input Parameters:
wfrPathName = '/home/andrei/Documents/9_term/diplom/BEAMLINE/files/' #example data sub-folder name
wfr1FileName = 'wfr1_for_BINP_UND.wfr' #file name for output UR flux data
wfr2FileName = 'wfr2_for_BINP_UND.wfr'
wfr3FileName = 'wfr3_for_BINP_UND.wfr'

stkFPathName = '/home/andrei/Documents/9_term/diplom/BEAMLINE/files/'
stkFFileName = 'stkF_for_BINP_UND.stkF'
#***********Undulator
undarr = []
undarrH = []
distz =  []
distx =  []
disty =  []
'''
Length = 2.3 #m
PER = 20
NumPIECE = 4
undper = 0.018
numper = PER/NumPIECE
magf = 1
magf_step = 10/100

Bx = 1 #Peak Horizontal field [T]
By = 1 + magf_step/2 #+ (NumPIECE-1)*magf_step#Peak Vertical field [T]
phBx = 0#Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx = 1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy = 1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 0 #Longitudinal Coordinate of Undulator Center [m]

und= []

for i in range(NumPIECE):
    B1 = Bx + i*magf_step
    und = SRWLMagFldU([SRWLMagFldH(1, 'v', B1, phBy, sBy, 1)], undper, numper)#, SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
    undarr.append(und)
    distz.append(-(NumPIECE - i - 1)*undper*numper)
    distx.append(0)
    disty.append(0)

for i in range(NumPIECE):
    B2 = By + i*magf_step
    und = SRWLMagFldU([SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undper, numper)#, SRWLMagFldH(1, 'h', B2, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
    undarr.append(und)
    distz.append((i+1)*undper*numper )
    distx.append(0)
    disty.append(0)
print(distx, distz)

'''
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
print("K = ",K)
#und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1), SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
#magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements
#magFldCnt = SRWLMagFldC(undarr, array('d', distx), array('d', disty), array('d', distz)) #Container of all Field Elements

#magFldCnt = SRWLMagFldC(undarr, array('d', distx), array('d', disty), array('d', distz)) #Container of all Field Elements
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
stkF = SRWLStokes() #for spectral flux vs photon energy
stkF.allocate(100, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.zStart = 22. #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = 1000. #initial photon energy [eV]
stkF.mesh.eFin = 1400. #final photon energy [eV]
a = 0.002
stkF.mesh.xStart = -a #initial horizontal position [m]
stkF.mesh.xFin = a #final horizontal position [m]
stkF.mesh.yStart = -a #initial vertical position [m]
stkF.mesh.yFin = a #final vertical position [m]


wfr1 = SRWLWfr() #For spectrum vs photon energy

wfr1.allocate(200, 40, 40) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 6800. #Initial Photon Energy [eV]
wfr1.mesh.eFin = 7200#4300. #Final Photon Energy [eV]
wfr1.avgPhotEn= 7008#4205
a = 0.002
wfr1.mesh.xStart = -a #Initial Horizontal Position [m]
wfr1.mesh.xFin = a #Final Horizontal Position [m]
wfr1.mesh.yStart = -a #Initial Vertical Position [m]
wfr1.mesh.yFin = a #Final Vertical Position [m]
wfr1.partBeam = eBeam

wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy

wfr2.allocate(1, 151, 151) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = 7008-40#4205 #Initial Photon Energy [eV]
wfr2.mesh.eFin = wfr2.mesh.eFin #Final Photon Energy [eV]
a = 0.001
wfr2.mesh.xStart = -a #Initial Horizontal Position [m]
wfr2.mesh.xFin = a #Final Horizontal Position [m]
wfr2.mesh.yStart = -a #Initial Vertical Position [m]
wfr2.mesh.yFin = a #Final Vertical Position [m]
wfr2.partBeam = eBeam
'''
wfr3 = SRWLWfr() #For spectrum vs photon energy

wfr3.allocate(200, 40, 40) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
#wfr3.avgPhotEn = 1260
wfr3.unitElFld = 2
#wfr3.unitElFldAng = 1
wfr3.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
#wfr3.mesh.eStart = -5.*10.0e-15 #Initial Time [s] 
#wfr3.mesh.eFin = 5.*10.0e-15 #Final Time [s]
wfr3.mesh.eStart = 6800#4100. #Initial Photon Energy [eV]
wfr3.mesh.eFin = 7200#4300. #Final Photon Energy [eV]
a = 0.002
wfr3.mesh.xStart = -a #Initial Horizontal Position [m]
wfr3.mesh.xFin = a #Final Horizontal Position [m]
wfr3.mesh.yStart = -a #Initial Vertical Position [m]
wfr3.mesh.yFin = a #Final Vertical Position [m]
wfr3.partBeam = eBeam
'''
#sys.exit(0)
#**********************Calculation (SRWLIB function calls)

print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
print('done')

print('   Performing Electric Field (wavefront at fixed photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr2, 0, magFldCnt, arPrecPar)
print('done')
'''
print('   Performing Electric Field (wavefront at fixed photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr3, 0, magFldCnt, arPrecPar)
print('done')

print('   Performing Spectral Flux (Stokes parameters) calculation ... ', end='')
srwl.CalcStokesUR(stkF, eBeam, und1, arPrecF)
print('done')
'''
#saving Wave Front to a file
afile = open(wfrPathName + wfr1FileName, 'wb')
pickle.dump(wfr1, afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'wb')
pickle.dump(wfr2, afile)
afile.close()
'''
afile = open(wfrPathName + wfr3FileName, 'wb')
pickle.dump(wfr3, afile)
afile.close()

afile = open(stkFPathName + stkFFileName, 'wb')
pickle.dump(stkF, afile)
afile.close()

'''












