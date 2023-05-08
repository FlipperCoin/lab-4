import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

L = 16e-3
W = 10e-3
d= 1e-3

#%%

Up = np.array([])
Ip = np.array([])

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
Ip = np.array([])

reg = linregress(Ip, UH)
RH = reg.slope*d/B

plt.errorbar(Ip, UH, 0, 0, '.', label='data')
plt.plot(smooth_I, reg.slope*smooth_I+reg.intercept, '--', label='fit')
plt.grid()
plt.legend()
plt.show()
