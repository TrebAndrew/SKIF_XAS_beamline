#############################################################################
#Extract two electric fields from files. Plot intensity distribution and spectrum.
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
#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.4 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = 3./0.51099890221e-03 #relative energy 3 Gev??
sigEperE = 0.00089 #relative RMS energy spread
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

'''
wfr1 = SRWLWfr() #For spectrum vs photon energy

wfr1.allocate(250, 40, 40) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 1200. #Initial Photon Energy [eV]
wfr1.mesh.eFin = 3500. #Final Photon Energy [eV]
a = 0.002
wfr1.mesh.xStart = -a #Initial Horizontal Position [m]
wfr1.mesh.xFin = a #Final Horizontal Position [m]
wfr1.mesh.yStart = -a #Initial Vertical Position [m]
wfr1.mesh.yFin = a #Final Vertical Position [m]
wfr1.partBeam = eBeam

wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy

wfr2.allocate(1, 151, 151) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = 22. #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = 2650 #Initial Photon Energy [eV]
wfr2.mesh.eFin = wfr2.mesh.eFin #Final Photon Energy [eV]
wfr2.mesh.xStart = -0.006 #Initial Horizontal Position [m]
wfr2.mesh.xFin = 0.006 #Final Horizontal Position [m]
wfr2.mesh.yStart = -0.006 #Initial Vertical Position [m]
wfr2.mesh.yFin = 0.006 #Final Vertical Position [m]
wfr2.partBeam = eBeam
'''
#sys.exit(0)
#**********************Calculation (SRWLIB function calls)

afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr2 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr3FileName, 'rb')
wfr3 =  pickle.load(afile)
afile.close()

afile = open(stkFPathName + stkFFileName, 'rb')
stkF =  pickle.load(afile)
afile.close()

            ######### Spectrum #######
print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
print('done')

print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])
'''
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI3 = array('f', [0]*wfr3.mesh.ne)
srwl.CalcIntFromElecField(arI3, wfr3, 6, 2, 0, wfr3.mesh.eStart, wfr3.mesh.xStart, wfr3.mesh.yStart)
print('done')

print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
uti_plot1d(arI3, [wfr3.mesh.eStart, wfr3.mesh.eFin, wfr3.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
'''
'''
uti_plot1d(stkF.arS, [stkF.mesh.eStart, stkF.mesh.eFin, stkF.mesh.ne], ['Photon Energy [eV]', 'Flux [ph/s/.1%bw]', 'Flux through Finite Aperture'])
'''
#uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])

            ######### Intensity #######
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
print('done')

uti_plot2d(arI2, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])
'''
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr1.mesh.eStart, 0, 0)
print('done')

arI2x = array('f', [0]*wfr2.mesh.nx) #array to take 1D intensity data (vs X)
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2x, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(horizontal cut at x = 0)'])

arI2y = array('f', [0]*wfr2.mesh.ny) #array to take 1D intensity data (vs Y)
srwl.CalcIntFromElecField(arI2y, wfr2, 6, 0, 2, wfr2.mesh.eStart, 0, 0)
uti_plot1d(arI2y, [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Vertical Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV\n(vertical cut at y = 0)'])
  ''' 
uti_plot_show() #show all graphs (and block execution)
print('done')


#***************** Optical Elements and Propagation Parameters
#***********Optical Elements and Propagation Parameters
#Sequence of Optical Elements:
#   <1st Crystal of DCM>
#   <2nd Crystal of DCM>
#   <Drift>

#C(400) Crystal Constants:
dSpC400 = 3.1355713563754857#(Si(111))#0.89178 #Crystal reflecting planes d-spacing for C(400) crystal
#psi0rC400 = -0.17530e-04; psi0iC400 = 0.21089e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
psi0rC400 = -3.389135132e-5; psi0iC400 = 1.650465185e-6#-0.21732e-04; psi0iC400 = 0.28005e-07 #Real and imaginary parts of 0-th Fourier component of crystal polarizability
#psihrC400 = -0.45300e-05; psihiC400 = 0.20314E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihrC400 = -1.797536216e-5; psihiC400 = 1.1523113e-6#-0.54377e-05; psihiC400 = 0.25934E-07 #Real and imaginary parts of h-th Fourier component of crystal polarizability
psihbrC400 = psihrC400; psihbiC400 = psihiC400 #Real and imaginary parts of -h-th Fourier component of crystal polarizability
thickCryst = 10e-03 #0.5e-03 #Thickness of each crystal [m]
angAsCryst = 0 #Asymmetry angle of each crystal [rad]

#1st Crystal:
opCr1 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst)
#Find appropriate orientation of the 1st crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr1 = opCr1.find_orient(wfr1.avgPhotEn-12)
orientCr1 = orientDataCr1[0] #1st crystal orientation
tCr1 = orientCr1[0]; nCr1 = orientCr1[2] # Tangential and Normal vectors to crystal surface
#print('   1st crystal orientation:'); print('   t=', tCr1, 's=', orientCr1[1], 'n=', nCr1)
#Set crystal orientation:
opCr1.set_orient(nCr1[0], nCr1[1], nCr1[2], tCr1[0], tCr1[1])
orientOutFrCr1 = orientDataCr1[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr1 = orientOutFrCr1[0]; ryCr1 = orientOutFrCr1[1]; rzCr1 = orientOutFrCr1[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   1st crystal output beam frame:'); print('   ex=', rxCr1, 'ey=', ryCr1, 'ez=', rzCr1)
TrM = [rxCr1, ryCr1, rzCr1] #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
uti_math.matr_print(TrM)

#2nd Crystal:
opCr2 = SRWLOptCryst(_d_sp=dSpC400, _psi0r=psi0rC400,
                     _psi0i=psi0iC400, _psi_hr=psihrC400, _psi_hi=psihiC400, _psi_hbr=psihbrC400, _psi_hbi=psihbiC400,
                     _tc=thickCryst, _ang_as=angAsCryst)
#Find appropriate orientation of the 2nd crystal and the corresponding output beam frame (in the incident beam frame):
orientDataCr2 = opCr2.find_orient( wfr1.avgPhotEn-12, _ang_dif_pl=3.141593)#E = 4205 wfr1.avgPhotEn

orientCr2 = orientDataCr2[0] #2nd crystal orientation
tCr2 = orientCr2[0]; nCr2 = orientCr2[2] # Tangential and Normal vectors to crystal surface
#print('   2nd crystal orientation:'); print('   t=', tCr2, 's=', orientCr2[1], 'n=', nCr2)
#Set crystal orientation:
opCr2.set_orient(nCr2[0], nCr2[1], nCr2[2], tCr2[0], tCr2[1])
orientOutFrCr2 = orientDataCr2[1] #Orientation (coordinates of base vectors) of the output beam frame 
rxCr2 = orientOutFrCr2[0]; ryCr2 = orientOutFrCr2[1]; rzCr2 = orientOutFrCr2[2] #Horizontal, Vertical and Longitudinal base vectors of the output beam frame
#print('   2nd crystal output beam frame:'); print('   ex=', rxCr2, 'ey=', ryCr2, 'ez=', rzCr2)
TrM = uti_math.matr_prod(TrM, [rxCr2, ryCr2, rzCr2]) #Input/Output beam transformation matrix (for debugging)
#print('   Beam frame transformation matrix (from the begining of opt. scheme to output of current element):')
uti_math.matr_print(TrM)
#print('   After the two crystals of DCM, the transformation matrix should be close to the unit matrix.')

slit = 1.5 # mm
slit_befour_MCh = 1.5#mm
SLIT = SRWLOptA('r', 'a', 2*slit*1e-03, 2*slit*1e-03)
Drift_SLIT_MCh = SRWLOptD(5.) #Drift from SSA to Center of Vertically Focusing K-B Mirror (VKB)

SLIT_BEFOUR_MCh = SRWLOptA('r', 'a', 2*slit_befour_MCh*1e-03, 2*slit_befour_MCh*1e-03)
Drift_AFTER_MCh = SRWLOptD(5.)
#Wavefront Propagation Parameters:
#                           [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
ppSLIT =                    [ 0,  0,  1.,  0,  0, 1.5, 1.0, 1.1,  6.,  0,  0,   0]

ppDrift_SLIT_MCh =          [ 0,  0,  1.,  0,  0, 1.0, 4.0, 1.0,  4.0, 0,  0,   0]
prPar0 =                    [ 0,  0,  1.,  1,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]
prParPost =                 [ 0,  0,  1.,  0,  0, 1.0, 1.0, 1.0,  1.0, 0,  0,   0]

#ppSLIT_BEFOUR_MCh =         [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
#ppDrift_AFTER_MCh =         [ 0,  0, 1.,   1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
#prParPost =                 [ 0,  0,  1.,  0,  0,  0.6, 4.5, 0.6, 4.5, 0, 0, 0]
#ppSSA =                     [0, 0, 1. ,0, 0, 1.1, 1.,   1.1,  1.0,0, 0, 0]
#prParInit =                 [0, 0, 1., 1, 0, 1.1, 2.,   1.1,  2., 0, 0, 0]
#prPar0 =    [0, 0, 1., 1, 0, 1.1, 1.,   1.1,  1., 0, 0, 0]

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
#optBL = SRWLOptC([SLIT,   Drift_SLIT_MCh,   SLIT_BEFOUR_MCh,   opCr1,  opCr2, Drift_AFTER_MCh], 
#                 [ppSLIT, ppDrift_SLIT_MCh, ppSLIT_BEFOUR_MCh, prPar0, prPar0,ppDrift_AFTER_MCh])
#optBL = SRWLOptC([opCr1,  opCr2], 
#                 [prPar0, prPar0, prParPost])
optBL = SRWLOptC([Drift_SLIT_MCh, opCr1,  opCr2], 
                 [ppDrift_SLIT_MCh, prPar0, prPar0, prParPost])
#///////////WAVEFRONT PROPOGATION//////#
'''
afile = open(wfrPathName + wfr1FileName, 'rb')
wfr1 =  pickle.load(afile)
afile.close()

afile = open(wfrPathName + wfr2FileName, 'rb')
wfr2 =  pickle.load(afile)
afile.close()
'''
print('   Simulating Electric Field Wavefront Propagation ... ', end='')
t0 = time.time();
srwl.PropagElecField(wfr2, optBL)
srwl.PropagElecField(wfr1, optBL)

print('done; lasted', round(time.time() - t0), 's')

print('   Extracting Intensity from the Propagated Electric Field  ... ', end='')
arI1 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" 2D array to take intensity data
srwl.CalcIntFromElecField(arI1, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
uti_plot2d(arI1, [1000*wfr2.mesh.xStart, 1000*wfr2.mesh.xFin, wfr2.mesh.nx], [1000*wfr2.mesh.yStart, 1000*wfr2.mesh.yFin, wfr2.mesh.ny], ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity at ' + str(wfr2.mesh.eStart) + ' eV'])

print('   Extracting Intensity from calculated Electric Field(Spectral Flux) ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 2, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
uti_plot1d(arI1, [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne], ['Photon Energy [eV]', 'Intensity [ph/s/.1%bw/mm^2]', 'On-Axis Spectrum'])






