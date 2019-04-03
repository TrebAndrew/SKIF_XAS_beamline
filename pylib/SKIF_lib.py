'''	SKIF lib
	@TrebAndrei
'''


import numpy as np
import matplotlib.pyplot as plt
from srwlib import *
from uti_plot import *

speed_of_light = 299792458 #[m/s]
h_bar = 6.582119514e-16 #[eV*s]
gamma = 3./0.51099890221e-03
e_ = 1.60218e-19 #elementary charge


def calc_bandwidth(x, y):
    y = y / np.sum(y)
    sigma = np.sqrt(np.sum((x)**2 * y))
    return(sigma)
    
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
        arI = arI/E/1e-3 * ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    elif elec_fld_units=='W/mm^2/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 0, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3 * ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    elif elec_fld_units is 'ph/s/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI/E/1e-3 * ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    elif elec_fld_units is 'W/eV':
        srwl.CalcIntFromElecField(arI, wfr, 6, 2, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
        arI = arI*E*e_/E/1e-3 * ((wfr.mesh.eFin-wfr.mesh.eStart)/wfr.mesh.ne)
    
    return E, arI
    
    
    
def skf_plot(x, y, color='blue', elec_fld_units=None, grid=False, log_x=False):
        
    plt.plot(x, y, color=color)
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
        plt.xlim(0)#, np.max(x))

    #plt.savefig('/home/andrei/Documents/SKIF_XAS_beamline/TeXDoc/pic/' +'spec_1_1'+'.pdf')
    #plt.show()

