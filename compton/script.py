""" Theory
Q1 
Descloizite (#10 - 𝑃𝑏𝑍𝑛(𝑂𝐻)𝑉𝑂4)
Element Kα1 Kα2 Kβ1 Lα1 Lα2 Lβ1 Lβ2 Lγ1 Mα1
82 Pb 74,969.4 72,804.2 84,936 10,551.5 10,449.5 12,613.7 12,622.6 14,764.4 2,345.5 
30 Zn 8,638.86 8,615.78 9,572.0 1,011.7 1,011.7 1,034.7
23 V 4,952.20 4,944.64 5,427.29 511.3 511.3 519.2 
8 O 524.9 

=> 
Pb 10,551.5 10,449.5 12,613.7 12,622.6 14,764.4
Zn 8,638.86 8,615.78 9,572.0
V 4,952.20 4,944.64 5,427.29

Q2
42 Mo 17,479.34 17,374.3 19,608.3 2,293.16 2,289.85 2,394.81 2,518.3 2,623.5
lmbda_tag - lmbda = h/(m_e*c) * (1-np.cos(theta))
=> lmbda_max = lmbda + 2*h/(m_e*c)
=> e_max = h*c/lmbda_max = h*c/(lmbda + 2*h/(m_e*c))

Q3
lmbda_e = h/(m_e*c)

Q4
PMMA (Plexiglass) - C5O2H8
Oxygen 532 eV
Carbon 283.8 eV
Hydrogen 13.6 eV 

Q5
approx 574.5 ev, from ~17420 to ~16850
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h,c,m_e,eV, physical_constants
from scipy.signal import find_peaks
import pandas as pd
r_e = physical_constants['classical electron radius'][0]

#%% Theory

# Q2
e_theo = np.array([17479.34, 17374.3, 19608.3])
lmbda_theo = h*c/(e_theo*eV)
lmbda_max = lmbda_theo + 2*h/(m_e*c)
e_min = h*c/lmbda_max/eV
max_delta_e =e_theo-e_min

# Q3
lmbda_e = h/(m_e*c)

# Q6
lmbda = 7.1e-11 # TODO
theta = np.linspace(0, np.pi, 50)
lmbda_tag = lmbda + h/(m_e*c) * (1-np.cos(theta))
diff_cross_sect = (1/2) * r_e**2 * (lmbda/lmbda_tag)**2 * ((lmbda/lmbda_tag) + (lmbda_tag/lmbda) - np.sin(theta)**2)
plt.figure(dpi=120)
plt.plot(np.rad2deg(theta), diff_cross_sect)
plt.grid()
plt.xlabel(r'$\theta^{\circ}$')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}$')
plt.show()

#%%
def hist_peaks(file, height=40, prominence=30, verb=False):
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts #running average
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    Mes1 = pd.read_csv(file,sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    peaks, _ = find_peaks(Counts_smoothed, height=height, prominence=prominence)
    plt.figure(dpi=150)
    plt.plot(Channels ,Counts_smoothed, label='Spectrum')
    plt.plot(Channels[peaks] ,Counts_smoothed[peaks],"x",color='red', 
    label='Lines')
    plt.grid()
    plt.ylabel('Impulses')
    plt.xlabel('Channels')
    plt.legend()
    if verb:
        return peaks, Channels, Counts, Counts_smoothed
    return peaks

def lin(x,a,b):
    return a*x+b
def get_e(a,b):
    def e(x):
        return lin(x, a, b) 
    return e

def hist_peaks_energy(file, height=40, prominence=30, verb=False):
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts #running average
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    Mes1 = pd.read_csv(file,sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    Energies = e(Channels)
    peaks, _ = find_peaks(Counts_smoothed, height=height, prominence=prominence)
    plt.figure(dpi=150)
    plt.plot(Energies ,Counts_smoothed, label='Spectrum')
    plt.plot(Energies[peaks] ,Counts_smoothed[peaks],"x",color='red', 
    label='Lines')
    plt.grid()
    plt.ylabel('Impulses')
    plt.xlabel('Energy')
    plt.legend()
    if verb:
        return peaks, Channels, Counts, Counts_smoothed
    return peaks, Counts_smoothed[peaks]

#%%

a = 8.193472040746805 # pm 0.030828000015000242
b = 952.863450735158 # pm 38.44782111295924

e = get_e(a,b)

#%%