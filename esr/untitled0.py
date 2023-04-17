# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:38:44 2023

@author: student
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, mu_0
import pandas as pd

# 

v0 = 98.6 # MHz, 0%
v1 = 98.3 # MHz, 13%
g = 2.0036
mu = 0.927e-23 # J/T

c_ratio_res = 18 # 18%

v_res = v0 + (v1 - v0)/13 * c_ratio_res

v_res *= 1e6

B_res = v_res*h/(g*mu)

p = pd.read_csv("D:/itai/min_mod_1 0.csv", sep=',', header=1, usecols=[3,4,5])

t = p['Time (s)']
v1 = p['1 (VOLT)']
v2 = p['2 (VOLT)']

plt.plot(t, v1, label='v1')
plt.plot(t, v2, label='v2')
plt.legend()
plt.grid()
plt.show()

V_amp1 = np.max(v1)
V_amp2 = np.min(v1)

V_amp = (V_amp1 + -1*V_amp2)/2
R_pdf = 0.82 
R_calc = 0.774
R=R_calc

I = V_amp/R 

K1 = B_res / (mu_0*I)

#%%

p = pd.read_csv("D:/itai/dc_1 0.csv", sep=',', header=1, usecols=[3,4,5])

I = 0.5077 # A, multimeter
V = 11.85 # Volt, around (11.8-11.9)
# A = 0.495-0.52 from dc supply


#%%

p = pd.read_csv("D:/itai/dc_2 0.csv", sep=',', header=1, usecols=[3,4,5])

I = 0.5046 # A, multimeter
V = 11.8 # Volt, around (11.75-11.85)
# A = 0.495-0.52 from dc supply

K2 = B_res / (mu_0 * I)

#%%

p = pd.read_csv("D:/itai/xy_1 0.csv", sep=',', header=1, usecols=[3,4,5])

I = 0.5200 # A, pm 0.0001
# power supply 0.513-0.527 A
V = 12.15 # Volts, pm 0.01

K3 = B_res / (mu_0 * I)


