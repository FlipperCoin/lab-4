# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.constants import h, mu_0
#from uncertainties import ufloat

#%%
def ufloat(a,b):
    return a

# %%
v0 = ufloat(98.6e6,0.1e6) # MHz, 0%
v1 = ufloat(98.3e6,0.1e6) # MHz, 13%
g = ufloat(2.0036, 0.0002)
mu = 0.927e-23 # J/T

c_ratio_res = ufloat(18,0.5) # 18%

v_res = v0 + (v1 - v0)/13 * c_ratio_res

B_res = v_res*h/(g*mu)

R_pdf = ufloat(0.82, 0.05*0.82)
R_calc = ufloat(402e-3,30e-3)/ufloat(0.5192,0.0001)
#0.774
R=R_pdf

# %%
p = pd.read_csv("data/min_mod_1 0.csv", sep=',', header=1, usecols=[3,4,5])

t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend(loc='upper center', bbox_to_anchor=(0.75,1.05))
plt.grid()
plt.show()

V_amp1 = np.max(v1)
V_amp2 = np.min(v1)

V_amp = (V_amp1 + -1*V_amp2)/2

I = V_amp/R 

K1 = B_res / (mu_0*I)
K1_high = K1

# %%
p = pd.read_csv("data/min_mod_2 0.csv", sep=',', header=1, usecols=[3,4,5])

t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

V_amp1 = np.max(v1)
V_amp2 = np.min(v1)

V_amp = (V_amp1 + -1*V_amp2)/2

I = V_amp/R 

K1_low = B_res / (mu_0*I)

plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend(loc='upper center', bbox_to_anchor=(0.75,1.00))
plt.grid()
plt.show()

# %%

p = pd.read_csv("data/min_mod_3 0.csv", sep=',', header=1, usecols=[3,4,5])

t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

V_amp1 = np.max(v1)
V_amp2 = np.min(v1)

V_amp = (V_amp1 + -1*V_amp2)/2

I = V_amp/R 

K1_mid = B_res / (mu_0*I)

plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend(loc='upper center', bbox_to_anchor=(0.75,1.00))
plt.grid()
plt.show()

# %%

p = pd.read_csv("data/dc_1 0.csv", sep=',', header=1, usecols=[3,4,5])
t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

I = 0.5077 # A, multimeter
V = 11.85 # Volt, around (11.8-11.9)
# A = 0.495-0.52 from dc supply

K2_1 = B_res / (mu_0 * I)
K2 = K2_1
plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend()
plt.grid()
plt.show()

right = np.max(v1)
left = np.min(v1)

mid = (right + left) / 2

# %%
plt.errorbar(v1, v2, 0.001, 0.001, '.', label='')
plt.xlabel(r'Outer Coil Voltage $\left[V\right]$')
plt.ylabel(r'Absorption Signal $\left[V\right]$')
plt.axvline(mid, color='orange', linestyle='dashed')
plt.grid()
plt.show()

# %%
p = pd.read_csv("data/dc_2 0.csv", sep=',', header=1, usecols=[3,4,5])
t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

I = ufloat(0.5046, 0.0001) # A, multimeter
V = ufloat(11.8, 0.05) # Volt, around (11.75-11.85)
# A = 0.495-0.52 from dc supply

K2_2 = B_res / (mu_0 * I)

plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend()
plt.grid()
plt.show()

right = np.max(v1)
left = np.min(v1)

mid = (right + left) / 2


# %%
plt.errorbar(v1, v2, 0.001, 0.001, '.', label='')
plt.xlabel(r'Outer Coil Voltage $\left[V\right]$')
plt.ylabel(r'Absorption Signal $\left[V\right]$')
plt.axvline(mid, color='orange', linestyle='dashed')
plt.grid()
plt.show()

# %%

p = pd.read_csv("data/xy_1 0.csv", sep=',', header=1, usecols=[3,4,5])
t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

plt.errorbar(t, v1, 0.001, 0, '.', label='Outer Coil')
plt.errorbar(t, v2, 0.001, 0, '.', label='Absorption')
plt.xlabel(r'Time $\left[sec\right]$')
plt.ylabel(r'Signal $\left[V\right]$')
plt.legend(bbox_to_anchor=(0.71,0.9))
plt.grid()
plt.show()

I = ufloat(0.5200, 0.0001) # A, pm 0.0001
# power supply 0.513-0.527 A
V = ufloat(12.15, 0.01) # Volts, pm 0.01

K3 = B_res / (mu_0 * I)

right = np.max(v1)
left = np.min(v1)

mid = (right + left) / 2


# %%

plt.errorbar(v1, v2, 0.001, 0.001, '.', label='')
plt.axvline(mid, color='orange', linestyle='dashed')
plt.xlabel(r'Outer Coil Voltage $\left[V\right]$')
plt.ylabel(r'Absorption Signal $\left[V\right]$')
plt.grid()
plt.show()


