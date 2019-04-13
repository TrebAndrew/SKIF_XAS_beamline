'''	SKIF lib
	@TrebAndrei
'''


import numpy as np
import matplotlib.pyplot as plt
from srwlib import *
from uti_plot import *
import SKIF_lib as skf


speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03
e_ = 1.60218e-19 #elementary charge


def calc_bandwidth(wfr):
    A = 1e6/wfr.mesh.zStart

    arIx = array('f', [0]*wfr.mesh.nx) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    arIx = arIx / np.sum(arIx)
    sigma_x = np.sqrt(np.sum((x)**2 * arIx))
    
    arIy = array('f', [0]*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    arIy = arIy / np.sum(arIy)
    sigma_y = np.sqrt(np.sum((x)**2 * arIy))
    return(sigma_x, sigma_y)
    
def pycry_trans(crystal='diamond', Emin=None, Emax=None, ne=None):
    
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
    if Emax==np.max(X_) and Emin==np.min(X_):
        step = (Emax-Emin)/ne
        for i in range(len(x_)-1):
            while round(x_[i] + j*step) < round(x_[i+1]):
                x.append(x_[i] + j*step)
                y.append(y_[i])
                j = j + 1
            j = 0
    
    return(np.array(y))
    
def renorm_wfr(wfr, elec_fld_units=None):
    arI = array('f', [0]*wfr.mesh.ne)
    E = np.linspace(wfr.mesh.eStart, wfr.mesh.eFin, wfr.mesh.ne)
    
    if elec_fld_units is None:
        srwl.CalcIntFromElecField(arI, wfr, 6, 0, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
    elif elec_fld_units is 'ph/s/mm^2/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 0, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    elif elec_fld_units=='W/mm^2/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 0, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3#*((wfr.mesh.eFin - wfr.mesh.eStart) / wfr.mesh.ne)
    elif elec_fld_units is 'ph/s/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    elif elec_fld_units is 'W/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3 #* ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    
    return E, arI  
    
    
def skf_plot(x, y, color='blue', elec_fld_units=None, grid=False, log_x=False, linewidth=1, save_fig=False, file_name=None):
    path_name = '/home/andrei/Documents/9_term/diplom/beamlines/1_1/'
    plt.plot(x, y, color=color, linewidth=linewidth)
    plt.xlabel(r'$E, [эВ]$', fontsize=14, labelpad = 0.0)
    
    if elec_fld_units is None:
        plt.ylabel(r'$a.u.$', fontsize=14, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'ph/s/mm^2/0.1%bw':
        plt.ylabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=14, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'ph/s/0.1%bw':
        plt.ylabel(r'$I, [\gamma/с/0.1\%ПП]$', fontsize=14, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'W/mm^2/eV':
        plt.ylabel(r'$I, [Вт/мм^2/эВ]$', fontsize=14, labelpad = 0.0, rotation=90)
    
    elif elec_fld_units == 'W/eV':
        plt.ylabel(r'$I, [Вт/эВ]$', fontsize=14, labelpad = 0.0, rotation=90)  
    
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
    ax.yaxis.set_label_coords(-0.07, 0.8)
    
    #plt.xticks([1200, 1300, 1400, 1500, 1600, 1700, 1800],fontsize=12)
            #  [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$'])
    #plt.xticks(fontsize=12)
    #plt.yticks(fontsize=12)
    #plt.yticks([1e14,2e14,3e14,4e14,5e14,6e14,7e14], fontsize=12)
    if grid is True:
        plt.grid()
    if log_x is True:
        ax.set_yscale('log')
        plt.ylim(0.1, np.max(y))
    else:
        plt.ylim(0)#, np.max(y))
        plt.xlim(0, np.max(x))

    if save_fig is True:
        plt.savefig(path_name + file_name, dpi=350)#, bbox_inches='tight')
        
    plt.show()

def skf_subplot_XY(wfr, save_fig=False, figure_name=None, units='mm'):
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
        A = 1e6/wfr.mesh.zStart
        xy_unit=r'$[\mu rad]$'

    plt.subplot(221)  
    cmap_ph = plt.get_cmap('viridis')
    arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    
    for j in range(len(y)):
        for i in range(len(x)):
            z.extend([arI[j*len(x) + i]])
        Z.append(z)
        z = []

    plt.pcolormesh(x, y, Z, cmap=cmap_ph)  
    plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=14, labelpad = 0.0, rotation=90)
    plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=14, labelpad = 0.0)
    plt.title('Intensity at ' + str(wfr.mesh.eStart) + ' eV', fontsize=14)
    plt.xlim(A*wfr.mesh.xStart,  A*wfr.mesh.xFin)
    plt.ylim(A*wfr.mesh.yStart, A*wfr.mesh.yFin)

    
    plt.subplot(223)
    arIx = array('f', [0]*wfr.mesh.nx) #array to take 1D intensity data (vs X)
    srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0)
    x = np.linspace(A*wfr.mesh.xStart, A*wfr.mesh.xFin, wfr.mesh.nx)
    plt.plot(x, arIx, color='blue')
    plt.ylabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=14, labelpad = 0.0, rotation=90)
    plt.xlabel(r'$Horizontal Position$' + xy_unit, fontsize=14, labelpad = 0.0)
    plt.grid(True)
    plt.xlim(A*wfr.mesh.xStart,  A*wfr.mesh.xFin)
    plt.ylim(0)
    
    plt.subplot(222)
    arIy = array('f', [0]*wfr.mesh.ny) #array to take 1D intensity data (vs X)
    srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, wfr.mesh.eStart, 0, 0)
    y = np.linspace(A*wfr.mesh.yStart, A*wfr.mesh.yFin, wfr.mesh.ny)
    plt.plot(arIy,y, color='blue')
    plt.xlabel(r'$I, [\gamma/с/мм^2/0.1\%ПП]$', fontsize=14, labelpad = 0.0, rotation=0)
    plt.ylabel(r'$Vertical Position$' + xy_unit, fontsize=14, labelpad = 0.0, rotation=90)
    plt.grid(True)
    plt.ylim(A*wfr.mesh.yStart, A*wfr.mesh.yFin)
    plt.xlim(0)

    if save_fig is True:
        plt.savefig(path_name + figure_name, dpi=150)#, bbox_inches='tight')
        
def skf_plot_spec(wfr1, crystal=True):
    
    plt.figure(figsize=(1.5*10,1.5*3))
    arI1 = array('f', [0]*wfr1.mesh.ne)
    srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
    E1 = np.linspace(wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne)
    E, spec = skf.renorm_wfr(wfr1, elec_fld_units='W/mm^2/eV')
    skf.skf_plot(E, spec, elec_fld_units='W/mm^2/eV')

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
        plt.show()










