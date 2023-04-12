# -*- coding: utf-8 -*-
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
#from uncertainties import ufloat
#from uncertainties import unumpy as unp

#%% 
Heater_Current = 0.26 #A 
fh1 = pd.read_csv('data/FH1_260.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 


maxs = []

for n in np.arange(1,len(Va1)-1):
    if (I1[n+1] - I1[n]) <= 0:
        if (I1[n-1] - I1[n]) <= 0:
            maxs = maxs + [n]    

maxs = np.array(maxs)
maxs = maxs[12:] # determined visually

va_maxs = Va1[maxs]
n = np.arange(len(va_maxs))
reg = linregress(n, va_maxs)


plt.figure()
plt.plot(n, va_maxs, '.', label='data')
plt.plot(n, reg.slope*n+reg.intercept, label='linear fit')
plt.grid()
plt.legend()
plt.show()

hg_level = reg.slope # deltaE
delta = reg.slope*(-1)+reg.intercept #phi_a-phi_k


#%%
Heater_Current = 0.23 #A 
fh1 = pd.read_csv('data/FH1_230.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 


maxs = []

for n in np.arange(1,len(Va1)-1):
    if (I1[n+1] - I1[n]) <= 0:
        if (I1[n-1] - I1[n]) <= 0:
            maxs = maxs + [n]    

maxs = np.array(maxs)
maxs = maxs[26:] # determined visually

va_maxs = Va1[maxs]
n = np.arange(len(va_maxs))
reg = linregress(n, va_maxs)


plt.figure()
plt.plot(n, va_maxs, '.', label='data')
plt.plot(n, reg.slope*n+reg.intercept, label='linear fit')
plt.grid()
plt.legend()
plt.show()

hg_level = reg.slope # deltaE
delta = reg.slope*(-1)+reg.intercept #phi_a-phi_k


#%%
Heater_Current = 0.24 #A 
fh1 = pd.read_csv('data/FH1_240ma_6200mv.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 


maxs = []

for n in np.arange(1,len(Va1)-1):
    if (I1[n+1] - I1[n]) <= 0:
        if (I1[n-1] - I1[n]) <= 0:
            maxs = maxs + [n]    

maxs = np.array(maxs)
maxs = maxs[13:] # determined visually

va_maxs = Va1[maxs]
n = np.arange(len(va_maxs))
reg = linregress(n, va_maxs)


plt.figure()
plt.plot(n, va_maxs, '.', label='data')
plt.plot(n, reg.slope*n+reg.intercept, label='linear fit')
plt.grid()
plt.legend()
plt.show()

hg_level = reg.slope # deltaE
delta = reg.slope*(-1)+reg.intercept #phi_a-phi_k


#%%
Heater_Current = 0.25 #A 
fh1 = pd.read_csv('data/FH1_250ma_6600mv.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 


maxs = []

for n in np.arange(1,len(Va1)-1):
    if (I1[n+1] - I1[n]) <= 0:
        if (I1[n-1] - I1[n]) <= 0:
            maxs = maxs + [n]    

maxs = np.array(maxs)
maxs = maxs[26:] # determined visually

va_maxs = Va1[maxs]
n = np.arange(len(va_maxs))
reg = linregress(n, va_maxs)


plt.figure()
plt.plot(n, va_maxs, '.', label='data')
plt.plot(n, reg.slope*n+reg.intercept, label='linear fit')
plt.grid()
plt.legend()
plt.show()

hg_level = reg.slope # deltaE
delta = reg.slope*(-1)+reg.intercept #phi_a-phi_k


#%%
Heater_Current = 0.24 #A 
fh1 = pd.read_csv('FH_Ion_Vr1_240ma_6200mv_final14v_step05_T110.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

c = Va1 > 9
Va1 = Va1[c]
I1 = I1[c]

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 

#%%
Heater_Current = 0.24 #A 
fh1 = pd.read_csv('FH_Ion_Vr1_240ma_6200mv_final15500mv_step05_T110.csv',sep='\t',header=5) # read the data. 
Va1 = np.array(fh1['Va(V)_1']) # accelerating voltage array 
I1 = np.array(fh1['Ia(E-12 A)_1']) # Current array 
T1 = np.array(fh1['T(c)_1']) #temperature array 

c = Va1 > 9
Va1 = Va1[c]
I1 = I1[c]

plt.figure() 
plt.plot(Va1, I1, label='Heater current {:.2f} [A]'.format(Heater_Current)) 
plt.ylabel('Current [pA]') 
plt.xlabel('Acceleration voltage [V]') 
plt.grid() 
