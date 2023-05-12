import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.constants import elementary_charge as e

L = 16e-3
W = 10e-3
d= 1e-3

#%%

Up = np.array([-1.668, -1.420, -1.127, -0.901, -0.628, -0.393, -0.133, 0.217, 0.433, 0.739,  1.019, 1.233, 1.495])
Ip = np.array([-30,    -25,    -20,    -15,    -10,    -5,     0,      5,     10,    15,     20,    25,    30])

reg = linregress(Ip, Up)
R0 = reg.slope

plt.errorbar(Ip, Up, 0, 0, '.', label='data')
smooth_I = np.linspace(-30e-3,30e-3,100)
plt.plot(smooth_I, R0*smooth_I+reg.intercept, '--', label='fit')
plt.grid()
plt.legend()
plt.show()

rhoxx = R0 / (L / (W*d))

#%%

B = 250e-3

UH = np.array([])
Ip = np.array([ ])

reg = linregress(Ip, UH)
RH = reg.slope*d/B

n = 1/(RH*e)

plt.errorbar(Ip, UH, 0, 0, '.', label='data')
plt.plot(smooth_I, reg.slope*smooth_I+reg.intercept, '--', label='fit')
plt.grid()
plt.legend()
plt.show()

#%%
