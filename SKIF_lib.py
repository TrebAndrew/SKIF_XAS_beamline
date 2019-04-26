'''	
SKIF lib
@author: Andrei Trebushini
'''

'''error debug functions must be added'''

import numpy as np
import matplotlib.pyplot as plt
from srwlib import *
from uti_plot import *
#import SKIF_lib as skf
import os

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03
e_ = 1.60218e-19 #elementary charge

def get_SKIF_directory():
    '''
    Function to get the project directory. Checked on Linux. 
    Probably does not work at the Window machines
    '''
    path = os.path.abspath(__file__)
    if path.endswith('SKIF_lib.py'):
        path = path[:-12]
    print('Your project name is ' + path + '\n', 'edit the skf.get_SKIF_directory if it is needed')
    return path

def calc_bandwidth(wfr, units):
    '''
    Calculate bandwidth of a given wavefront
    :wfr:    SRWLWfr object 
    :units: units of x and y coordinates mm or urad
    ''' 
    if units == 'mm': #reperesentation in mm
        A = 1e3
        xy_unit=r'$[mm]$'
    elif units == 'urad':#reperesentation in angles
        A = 1e6/wfr.mesh.zStart
        xy_unit=r'$[\mu rad]$'
    
    arIx = array('f', [0]*wfr.mesh.nx) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIx, wfr, 6, 1, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    arIx = arIx / np.sum(arIx) #normalise the array on unity
    sigma_x = np.sqrt(np.sum((x)**2 * arIx)) #std formula
    
    arIy = array('f', [0]*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIy, wfr, 6, 1, 2, wfr.mesh.eStart, 0, 0)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    arIy = arIy / np.sum(arIy) #normalise the array on unity
    sigma_y = np.sqrt(np.sum((x)**2 * arIy)) #std formula
    
    return(sigma_x, sigma_y)
    
def calc_FWHM(wfr, units):
    '''
    Calculate bandwidth of a given wavefront
    :wfr:    SRWLWfr object 
    :units: units of x and y coordinates mm or urad
    ''' 
    if units == 'mm': 
        A = 1e3
        xy_unit=r'$[mm]$'
    elif units == 'urad': #reperesentation in angles
        A = 1e6/wfr.mesh.zStart
        xy_unit=r'$[\mu rad]$'
    
    arIx = array('f', [0]*wfr.mesh.nx) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIx, wfr, 6, 1, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    difference = np.max(arIx) - np.min(arIx)
    HM = difference / 2
    nearest = (np.abs(arIx - HM)).argmin()
    max_index = np.max(np.where(arIx == np.amax(arIx)))
    FWHM_x = 2*(np.abs(x[nearest] - x[max_index]))
    
    arIy = array('f', [0]*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIy, wfr, 6, 1, 2, wfr.mesh.eStart, 0, 0)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    difference = np.max(arIy) - np.min(arIy)
    HM = difference / 2
    nearest = (np.abs(arIy - HM)).argmin()
    max_index = np.max(np.where(arIy == np.amax(arIy)))
    FWHM_y = 2*(np.abs(y[nearest] - y[max_index]))
    
    return(FWHM_x, FWHM_y)

def find_sigma_for_array(Y, X):
    '''
    Calculate bandwidth of a given array
    :Y, X:    arrays to find fwhm
    '''
    Y = Y / np.sum(Y) #normalise the array on unity
    X_mean = np.sum(X*Y)
    sigma = np.sqrt(np.sum((X-X_mean)**2 * Y)) #std formula
    return sigma

def find_FWHM_for_array(Y, X):
    '''
    Calculate bandwidth of a given array
    :Y, X:    arrays to find fwhm
    '''
    difference = np.max(Y) - np.min(Y)
    HM = difference / 2
    nearest = (np.abs(Y - HM)).argmin()
    max_index = np.max(np.where(Y == np.amax(Y)))
    FWHM = (np.abs(X[nearest] - X[max_index]))
    return FWHM

def pycry_trans(crystal='diamond', Emin=None, Emax=None, ne=None):
    '''
    Extract data for crystals absorption curves from your local repository. In development...
    :crystal: crystal data to be extracted
    :Emin:    lift bound of the energy 
    :Emax:    right bound of the energy 
    :ne:      number of points for the energy range
    '''
    pathname = '/home/andrei/Documents/SKIF_XAS_beamline/crystals_data/diamond_T/'

    if crystal is 'diamond':
        filename_x = 'diamond_T_x.txt'
        filename_y = 'diamond_T_y.txt'
        
    file_x = pathname + filename_x
    file_y = pathname + filename_y
    f = open(file_x,"r+")
    lines=f.readlines()
    X=[]
    
    for x in lines:
        X.append(x.split('\t')[0])
    f.close()
    
    f = open(file_y,"r+")
    lines=f.readlines()
    Y=[]
    
    for x in lines:
        Y.append(x.split('\t')[0])
    f.close()
        
    X_ = [0]*(len(X)-1)
    Y_ = [0]*(len(X)-1)
    for i in range(len(X)-1):
        X_[i]=float(X[i])
        Y_[i]=float(Y[i])
    
    x_ = []
    x_ = X_
    y_ = []
    y_ = Y_
    j = 0
    x = []
    y = []
    if Emax==np.max(X_) and Emin==np.min(X_) and ne != len(X):
        step = (Emax-Emin)/ne
        for i in range(len(x_)-1):
            while round(x_[i] + j*step) < round(x_[i+1]):
                x.append(x_[i] + j*step)
                y.append(y_[i])
                j = j + 1
            j = 0
    #elif Emax==np.max(X_) and Emin==np.min(X_) and ne==len(X):
        
    return(np.array(y))

def read_crystal_trans(file_path=None, file_name=None):
    
    if file_path is None:
        file_path = '/home/andrei/Documents/SKIF_XAS_beamline/crystals_data/diamond_T/diamond_T_100.txt'

    f = open(file_path+file_name,"r")
    lines=f.readlines()
    X=[]
    Y=[]
    for x in lines:
        X.append(float(x.split('\t')[0]))
        Y.append(float(x.split('\t')[1]))
    f.close()
    
    return(np.array(X), np.array(Y))


def renorm_wfr(wfr, elec_fld_units=None, emittance=0):
    '''
    Renormalise the wavefront from ph/s/mm^2/0.1%bw or ph/s/0.1%bw on desired units ph/s/mm^2/eV, W/mm^2/eV or
    or integrated over some aperture ph/s/eV, W/eV
    :wfr:    SRWLWfr object 
    :elec_fld_units: units to be transformed from ph/s/mm^2/eV, W/mm^2/eV, ph/s/eV, W/eV
    :emittance: 0 - filament beam, 1 - finite emittance beam
    '''
    arI = array('f', [0]*wfr.mesh.ne)
    E = np.linspace(wfr.mesh.eStart, wfr.mesh.eFin, wfr.mesh.ne)
    
    if elec_fld_units is None:
        srwl.CalcIntFromElecField(arI, wfr, 6, emittance, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)

    elif elec_fld_units is 'ph/s/mm^2/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, emittance, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    
    elif elec_fld_units=='W/mm^2/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, emittance, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3#*((wfr.mesh.eFin - wfr.mesh.eStart) / wfr.mesh.ne)
    
    elif elec_fld_units is 'ph/s/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2 + emittance, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    
    elif elec_fld_units is 'W/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2 + emittance, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    
    return E, arI  
    
    
def skf_plot(x, y, color='blue', scale_x=None, elec_fld_units=None, grid=True, log_x=False, 
             linewidth=1, save_fig=False, figure_name=None, show=True, file_path=None):
    '''
    Plot a graph with x and y. Optimised for spectra plotting although may be used for other dependencies.
    :x: 1d array 
    :y: 1d array
    :color: color of the line
    :elec_fld_units: SRWLWfr object units
    :grid: grid on/off
    :log_x: log scale on/off
    :linewidth: width of the line
    :save_fig:  flag for saving the graph
    :figure_name: filename of the saved file
    :show:      draw the graph of not
    '''
    #plt.figure(figsize=(1.5*4,1.5*3))
    if file_path is None:
        file_path = '/home/andrei/Documents/9_term/diplom/beamlines/1_4/'
    if scale_x is not None:
        x = x/scale_x
        plt.xlabel(r'$E, [кэВ]$', fontsize=18, labelpad = 0.0)
    else:
        plt.xlabel(r'$E, [эВ]$', fontsize=18, labelpad = 0.0)
    plt.plot(x, y, color=color, linewidth=linewidth)
    
    
    if elec_fld_units is None:
        plt.ylabel(r'$a.u.$', fontsize=18, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'ph/s/mm^2/0.1%bw':
        plt.ylabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=18, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'ph/s/0.1%bw':
        plt.ylabel(r'$I, [\gamma/с/0.1\%ПП]$', fontsize=18, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'W/mm^2/eV':
        plt.ylabel(r'$I, [Вт/мм^2/эВ]$', fontsize=18, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'W/eV':
        plt.ylabel(r'$I, [Вт/эВ]$', fontsize=18, labelpad = 0.0, rotation=90)  

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_position(('axes',0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_position(('data',0))
    ax.xaxis.set_label_coords(0.95, -0.12)
    ax.yaxis.set_label_coords(-0.1, 0.8)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    if grid is True:
        plt.grid()
    
    if log_x is True:
        ax.set_yscale('log')
        plt.ylim(0.1, np.max(y))
    else:
        plt.xlim(np.min(x), np.max(x))
        plt.ylim(0, np.max(y))
    plt.tight_layout()
    if save_fig is True:
        plt.savefig(file_path + figure_name)#, bbox_inches='tight')
        
#    if show is True:
#        plt.show()

def skf_wfr_subplot_XY(wfr, save_fig=False, figure_name=None, units='urad', fourth_plot=None, three_first=1, show=True, file_path=None):
    '''
    In development
    Draw a plot with four subplots with wfr field distribution. The fourth is for arbitrary dependence.
    :wfr:    SRWLWfr object 
    :save_fig:  flag for saving the graph
    :figure_name: filename of the saved file
    :units: units of x and y coordinates mm or urad
    :fourth_plot: flag for electric field calculation for four subplot
 *             0- "Single-Electron" Intensity; 
 *             1- "Multi-Electron" Intensity; 
 *             2- "Single-Electron" Flux; 
 *             3- "Multi-Electron" Flux; 
 *             4- "Single-Electron" Radiation Phase; 
 *             5- Re(E): Real part of Single-Electron Electric Field;
 *             6- Im(E): Imaginary part of Single-Electron Electric Field;
 *             7- "Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence);
    :three_first: flag for electric field calculation for four subplot (see above) 
    !!!be careful not all flags work properly and axis labels must be modified at the appropriate ones!!!
    :show:      draw the graph of not
    '''
    if file_path is None:
        file_path = '/home/andrei/Documents/9_term/diplom/beamlines/1_1/'
    z = []
    Z = []
    i = 0 
    j = 0   
    plt.figure(figsize=(1.5*8,1.5*8))
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    
    if units == 'mm': 
        A = 1e3
        xy_unit=r'$[mm]$'
    elif units == 'urad':
        A = 1e6/wfr.mesh.zStart
        xy_unit=r'$[\mu rad]$'
   
    sigma_x_mm, sigma_y_mm   = skf.calc_bandwidth(wfr, units='mm')
    sigma_x_rad, sigma_y_rad = skf.calc_bandwidth(wfr, units='urad')
#    print('sigma_x = ', round(sigma_x_rad,3),'[urad] \t','sigma_y = ', round(sigma_y_rad,3),'[urad]')
#    print('sigma_x = ', round(sigma_x_mm,3),'[mm] \t','sigma_y = ', round(sigma_y_mm,3),'[urad]')

    ########
    plt.subplot(221)  
    cmap_ph = plt.get_cmap('viridis')
    arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, three_first, 3, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    
    for j in range(len(y)):
        for i in range(len(x)):
            z.extend([arI[j*len(x) + i]])
        Z.append(z)
        z = []
    
    plt.pcolormesh(x, y, Z, cmap=cmap_ph)  
    plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=18, labelpad = 0.0, rotation=90)
    plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=18, labelpad = 0.0)
    plt.title('Multiple electron Intensity at ' + str(wfr.mesh.eStart) + ' eV', fontsize=18)
    plt.xlim(A*wfr.mesh.xStart,  A*wfr.mesh.xFin)
    plt.ylim(A*wfr.mesh.yStart, A*wfr.mesh.yFin)
    
    ########
    plt.subplot(224) 

    plt.axis('off')
    fontsize = 24
    x_pos = 0.1#*A*wfr.mesh.xFin
    y_pos = 0.9#*np.max(arIx)
    a = plt.text(x_pos, y_pos, r'$RMS_x = {} [\mu rad]$'.format(round(sigma_x_rad, 1)), fontsize=fontsize)
    y_pos = 0.75#*np.max(arIx)
    a = plt.text(x_pos, y_pos, r'$RMS_x = {} [mm]$'.format(round(sigma_x_mm, 3)), fontsize=fontsize)   
    y_pos = 0.6#*A*wfr.mesh.yFin
    b = plt.text(x_pos, y_pos, r'$RMS_y = {} [\mu rad]$'.format(round(sigma_y_rad, 1)), fontsize=fontsize)
    y_pos = 0.45#*A*wfr.mesh.yFin
    b = plt.text(x_pos, y_pos, r'$RMS_y = {} [mm]$'.format(round(sigma_y_mm, 3)), fontsize=fontsize)

    ########
    plt.subplot(223)
    arIx = array('f', [0]*wfr.mesh.nx) #array to take 1D intensity data (vs X)
    srwl.CalcIntFromElecField(arIx, wfr, 6, three_first, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    plt.plot(x, arIx, color='blue')
    plt.ylabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=18, labelpad = 0.0, rotation=90)
    plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=18, labelpad = 0.0)
    plt.grid(True)
    plt.xlim(A*wfr.mesh.xStart,  A*wfr.mesh.xFin)
    plt.ylim(0)

    ########
    plt.subplot(222)
    arIy = array('f', [0]*wfr.mesh.ny) #array to take 1D intensity data (vs X)
    srwl.CalcIntFromElecField(arIy, wfr, 6, three_first, 2, wfr.mesh.eStart, 0, 0)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    plt.plot(arIy, y, color='blue')
    plt.xlabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=18, labelpad = 0.0, rotation=0)
    plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=18, labelpad = 0.0, rotation=90)
    plt.grid(True)
    plt.ylim(A*wfr.mesh.yStart, A*wfr.mesh.yFin)
    plt.xlim(0)
    
    plt.tight_layout()
    if save_fig is True:
        plt.savefig(file_path + figure_name, dpi=200)#, bbox_inches='tight')
    if show is True:
        plt.show()
    else: 
        return None
        
def skf_power_subplot_XY(stkP, wfr=None,save_fig=False, figure_name=None, units='urad', show=True, path_name=None):
    '''
    In development
    Draw a plot with four subplots with power density distribution. The fourth is for arbitrary dependence (not implemented).
    :stkP:    SRWLStokes object 
    :save_fig:  flag for saving the graph
    :figure_name: filename of the saved file
    :show:      draw the graph of not
    '''
    if path_name is None:
        path_name = '/home/andrei/Documents/9_term/diplom/beamlines/1_1/'
    z = []
    Z = []
    i = 0 
    j = 0 

    plt.figure(figsize=(1.5*8,1.5*8))
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    if units == 'mm': 
        A = 1e3
        xy_unit=r'$[mm]$'
    elif units == 'urad':
        A = 1e6/stkP.mesh.zStart
        xy_unit=r'$[\mu rad]$'
    else:
        print('\n !Error! units must be mm or urad')
        return 0
    
    ########
    plt.subplot(221)  
    cmap_ph = plt.get_cmap('viridis')
    
    x = np.linspace(A*stkP.mesh.xStart, A*stkP.mesh.xFin, stkP.mesh.nx)
    y = np.linspace(A*stkP.mesh.yStart, A*stkP.mesh.yFin, stkP.mesh.ny)
    
    for j in range(len(y)):
        for i in range(len(x)):
            z.extend([stkP.arS[j*len(x) + i]])
        Z.append(z)
        z = []
    Z = np.array(Z)
    plt.pcolormesh(x, y, Z, cmap=cmap_ph)  
    plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=18, labelpad = 0.0, rotation=90)
    #plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=18, labelpad = 0.0)
    plt.title('Power density', fontsize=18)
    plt.ylim(A*stkP.mesh.yStart, A*stkP.mesh.yFin)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    ########    
    plt.subplot(223)
    x = np.linspace(A*stkP.mesh.xStart, A*stkP.mesh.xFin, stkP.mesh.nx)
    plt.plot(x, Z[int(stkP.mesh.nx*0.5),:], color='blue')
    plt.ylabel(r'$I, [Вт/мм^2]$', fontsize=18, labelpad = 0.0, rotation=90)
    plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=18, labelpad = 0.0)
    plt.grid(True)
    plt.xlim(A*stkP.mesh.xStart,  A*stkP.mesh.xFin)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim(0)

    ########    
    plt.subplot(222)
    y = np.linspace(A*stkP.mesh.yStart, A*stkP.mesh.yFin, stkP.mesh.ny)
    plt.plot(Z[:,int(stkP.mesh.ny*0.5)], y, color='blue')
    plt.xlabel(r'$I, [Вт/мм^2]$', fontsize=18, labelpad = 0.0, rotation=0)
    #plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=18, labelpad = 0.0, rotation=90)
    plt.grid(True)
    plt.ylim(A*stkP.mesh.yStart, A*stkP.mesh.yFin)
    plt.xlim(0)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    a = plt.text(np.max(Z)*0.15, A*stkP.mesh.yFin*0.8, r'$W_m = {} [W/mm^2]$'.format(round(np.max(stkP.arS),1)), fontsize=24)
    #########
    if wfr is not None:
        plt.subplot(224) 
        E, spec = skf.renorm_wfr(wfr, elec_fld_units=None, emittance=0)
        skf.skf_plot(E, spec, elec_fld_units='ph/s/mm^2/0.1%bw', scale_x = 1000)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)  
    print('max power density on axis = ', round(np.max(stkP.arS),1), '[W/mm^2]')

    if save_fig is True:
        plt.savefig(path_name + figure_name, dpi=150)#, bbox_inches='tight')
    
    if show is True:
        plt.show()
   
def skf_plot_spec(wfr1, dep_type=0, crystal=False):
    
    plt.figure(figsize=(1.5*10,1.5*3))
    E, spec = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV', emittance=dep_type)
    skf.skf_plot(E, spec, elec_fld_units='W/mm^2/eV')
    plt.show()
    if crystal is True:
        T = skf.pycry_trans(crystal='diamond', Emin=wfr1.mesh.eStart, Emax=wfr1.mesh.eFin, ne=wfr1.mesh.ne)
        spec = spec[:len(T)]
        E = E[:len(T)]
        skf.skf_plot(E, spec*T, color='red', elec_fld_units='W/mm^2/eV', grid=True)
        
        W = np.sum(spec) - np.sum(spec*T)

        plt.text(0.7, 0.89, r'$W_a = {} \; [Вт/мм^2]$'.format(round(W),1), {'color': 'black', 'fontsize': 20},
        horizontalalignment='left',
        verticalalignment='center',
        rotation=0,
        clip_on=False,
        transform=plt.gca().transAxes)
#       plt.show()











