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
approx 574.5 ev, from ~17426 to ~16850
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h,c,m_e,eV, physical_constants
from scipy.signal import find_peaks
from scipy.stats import linregress
import pandas as pd
from uncertainties import ufloat
from uncertainties import unumpy as unp
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
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

chnls = hist_peaks('data/Mo', height=10, prominence=10)


#%%

chnls = hist_peaks('data/Descloizite', height=5, prominence=10)

#%%
chan=np.array([2027, 2294, 955, 1179, 1435])
energies=np.array([17420, 19608, (8638.86 + 8615.78)/2, (10551.5 + 10449.5)/2, (12613.7 + 12622.6)/2])

reg=linregress(chan,energies)
xchan=np.arange(750,2500)
plt.figure(dpi=300)
plt.errorbar(chan ,energies,0,0,fmt='.', label='Lines')
plt.plot(xchan, reg.slope*xchan + reg.intercept, '--', label='Fit')
plt.grid()
plt.ylabel('Energy [eV]')
plt.xlabel('Channels')
plt.legend()
plt.show()

a = ufloat(reg.slope, reg.stderr)
b = ufloat(reg.intercept, reg.intercept_stderr)

#%%
def conv(channel):
    return a*channel + b

#%%

import pandas as pd # read the data from files
import matplotlib.pyplot as plt # plot the data
import numpy as np # just for fun :)
from scipy.signal import find_peaks # find the first estimate to a line.
from scipy.stats import linregress # for calibration and compton fit
from scipy.optimize import curve_fit # for gausian fit 
import scipy.constants as const # physical constants. 

#Define the Gaussian function
def Gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

#Define the Gaussian function
def uGauss(x, H, A, x0, sigma):
    return H + A * unp.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

# Mean over near chanels.
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

# just linear function
def lin(x,a,b):
    return a*x+b

# convertion of energies in eV to wavelength in m
def e2lam(e):
    hc=1.24*10**-6; # eV*m
    return hc/e

# compton shift in whavelength
def comp(theo_angles):
    h=6.626*10**-34 # eV*s
    c=2.9*10**8; # m/s
    me=9.1*10**-31; # kg
    return h/(me*c)*(1-np.cos(theo_angles))+e2lam(17426)

# The Klein-Nishina formula to energy 17426 eV
def klein_nish(A,theo_angles):
    E0=17426 # eV . The line energy without the Plexiglas
    lam0=e2lam(E0)
    return 1/2*A**2*(lam0/comp(theo_angles))**2*(lam0/comp(theo_angles)+comp(theo_angles)/lam0-np.sin(theo_angles)**2)

# Define the line energy and amplitudes
comp_amp=[]
comp_eng=[]
# Define 14 colors
colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 
'pink', 'brown', 'black', 'gray', 'cyan', 'magenta', 'lime', 'navy']
# The loop run on different angles
plt.figure(dpi=300)
for i in range(1,15):
    # import the data 
    data=pd.read_csv('data/run{}'.format(i), sep='\t', header=1)
    chanals=np.array(data['Channel/#'])
    Impulses=np.array(data['Impulses/#'])
    Imp_smooth=smooth(Impulses,50)
    
    # cut relevant interval
    x=noms(conv(chanals))
    y=Imp_smooth
    indss = (x>16000) & (x<18500)
    x=x[indss]
    y=y[indss]
    
    # plot the relevant interval
    plt.plot(x,y,':',color=colors[i-1], label=f'Run {i}')
    
    # first estimate the line energy
    peaks, properties = find_peaks(y, prominence=5, width=20,distance=1000)
    plt.plot(x[peaks], y[peaks], "rx") # plot the estimation 
 
    # Fit the line to gaussian. p0 is the initial guess 
    parameters, covariance = curve_fit(Gauss, x, y,p0=[10,y[peaks].item(),x[peaks].item(), 500])
    parameters_errs = np.sqrt(np.diag(covariance))
    uparams = unp.uarray(parameters, parameters_errs)

    plt.plot(x,noms(uGauss(x,uparams[0],uparams[1],uparams[2],uparams[3])),'--',color=colors[i-1])
    #acumulate the line energies and amplitudes
    comp_amp.append(uGauss(uparams[2],uparams[0],uparams[1],uparams[2],uparams[3]))
    comp_eng.append(uparams[2]) # e

plt.xlabel(r'Energy [eV]')
plt.ylabel(r'Intensity [au]')
plt.grid()
plt.legend()
plt.show()

#%%

comp_amp=np.array(comp_amp)
comp_eng=np.array(comp_eng)
angles = np.arange(20, 150+1, 10)
lmbda = e2lam(comp_eng)

x =(1-np.cos(np.deg2rad(angles)))[:len(lmbda)]
y = lmbda
reg = linregress(x, noms(y))

plt.figure(dpi=120)
plt.errorbar(x, noms(y), devs(y), 0, fmt='.', label='data')
plt.plot(x, reg.slope*x+reg.intercept, '--', label='fit')
plt.xlabel(r'$1-\cos\left(\theta\right)$')
plt.ylabel(r"$\lambda'$ $\left[m\right]$")
plt.legend()
plt.grid()
plt.show()

h=6.626*10**-34 # eV*s
c=2.9*10**8; # m/s
m_e = 1/ufloat(reg.slope, reg.stderr) * (h/c)

#%%
smooth_angles = np.linspace(0, 180, 100)
A=11
plt.figure(dpi=120)
plt.plot(smooth_angles, klein_nish(A, np.deg2rad(smooth_angles))-32, label='Klein-Nishina')
plt.errorbar(angles, noms(comp_amp), devs(comp_amp), 0, fmt='.', label='Data')
plt.xlabel(r'$\theta$ $\left[Degree\right]$')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}$')
plt.legend()
plt.grid()
plt.show()

#%%
def fit_klein_nish(theta, A):
    return klein_nish(A, theta)

res = curve_fit(fit_klein_nish, np.deg2rad(noms(angles)), noms(comp_amp))
print(f"A: {ufloat(res[0],np.sqrt(res[1][0]))}")
plt.plot(smooth_angles, klein_nish(res[0], np.deg2rad(smooth_angles)), '--')
plt.plot(noms(angles), noms(comp_amp), '.')
plt.grid()

# %%

chnls = hist_peaks('data/run5', height=5, prominence=10)

# %%

i=8
data=pd.read_csv('data/run{}'.format(i), sep='\t', header=1)
chanals=np.array(data['Channel/#'])
Impulses=np.array(data['Impulses/#'])
Imp_smooth=smooth(Impulses,10)

# cut relevant interval
x=noms(conv(chanals))
y=Imp_smooth
indss = (x>16000) & (x<18500)
x=x[indss]
y=y[indss]

plt.figure(dpi=150)
# plot the relevant interval
plt.plot(x,y,':', label=f'Run {i}')

# first estimate the line energy
peaks, properties = find_peaks(y, prominence=5, width=20,distance=1000)
plt.plot(x[peaks], y[peaks], "rx") # plot the estimation 

# Fit the line to gaussian. p0 is the initial guess 
parameters, covariance = curve_fit(Gauss, x, y,p0=[10,y[peaks].item(),x[peaks].item(), 500])
parameters_errs = np.sqrt(np.diag(covariance))
uparams = unp.uarray(parameters, parameters_errs)

# plt.plot(x,noms(uGauss(x,uparams[0],uparams[1],uparams[2],uparams[3])),'--')
print(uparams)