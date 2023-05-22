# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, poisson, norm
from uncertainties import ufloat
from uncertainties import unumpy as unp
from uncertainties.unumpy import uarray
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
from scipy.special import factorial

# %% [markdown]
# ### Statistics of counting and Background Radiation Measurement

# %%
# range 500 - 1200 v, steps 15v, 3 seconds per step, 
plat_df = pd.read_csv('data/plat.tsv', sep='\t+', header=9)

rate = plat_df['Counts'] / 3
plt.figure(dpi=100)
plt.errorbar(plat_df['Voltage'], rate, xerr=1, yerr=0, fmt='.')
plt.xlabel('Voltage [V]')
plt.ylabel(r'cps $\left[\frac{1}{s}\right]$')
plt.grid()
plt.show()

op_v = 950 # v

# %%
background_count = 43
background_rate = background_count / 100 # per sec

# %%
init_nbar = np.sum([0,5,5,5,2,3,3,1,4,2]) / 10
m = 150*init_nbar

# %%
df = pd.read_csv('data/nbar3.tsv', sep='\t+', header=9)

# %%
n = df['Counts']
n -= background_rate

# %%
print(f"m = {m} and len(n) = {len(n)}")

nbar = np.mean(n)
n_std = np.std(n)
nbar_std = n_std/np.sqrt(len(n)-1)
# n_std = np.sqrt(np.mean((n-nbar)**2)) # *np.sqrt(m)

K3 = (n-nbar)**3
K3_bar = np.mean(K3)
K3_std = np.std(K3)
K3_bar_std = K3_std/np.sqrt(len(n)-1)
# K3_std = np.sqrt(np.mean((K3-K3_bar)**2)) # *np.sqrt(m)

print(f"n +/- STD(nbar) = {nbar} +/- {nbar_std}")
print(f"K3 +/- STD(K3) = {K3_bar} +/- {K3_bar_std}")

bins=np.array([0,1,2,3,4,5,6,7,8])-background_rate
plt.figure(dpi=130)
plt.hist(n,bins)
plt.xlabel('counts in 1 sec intervals')
plt.ylabel('occurrences')
plt.grid()
plt.show()


# %% [markdown]
# ### Inverse Square Law

# %%
shelves_distance = ufloat(1.09e-2,0.01e-2) # guide ~10e-3
top_shelf_GM_window = ufloat(1.09e-2,0.01e-2) # guide 12.3e-3
# from second shelf 1.94e-3 
stage_thickness = 1e-3
number_of_shelves = 10

time_per_shelf = np.array([ 60,  50,    50,   40,    30,    20,    20,  8,    5,   3])
count_per_shelf = np.array([1320,1198 ,1546, 1527, 1493,  1312,   1936,  1154, 1209,1303])

rate_per_shelf = uarray(count_per_shelf,np.sqrt(count_per_shelf)) / time_per_shelf
normal_rate_per_shelf = rate_per_shelf - background_rate

distance_per_shelf = (number_of_shelves-1)*(shelves_distance) + top_shelf_GM_window - np.arange(number_of_shelves)*shelves_distance

plt.errorbar(noms(distance_per_shelf), noms(normal_rate_per_shelf), xerr=devs(distance_per_shelf),yerr=devs(normal_rate_per_shelf),fmt='.')
plt.grid()
plt.show()

x = distance_per_shelf
y = 1/unp.sqrt(normal_rate_per_shelf)
fit = linregress(noms(x), noms(y))

a = ufloat(fit.intercept,fit.intercept_stderr)/ufloat(fit.slope, fit.stderr)

# %%
smooth_x = np.linspace(np.min(noms(x)), np.max(noms(x)), 100)
plt.figure(dpi=150)
plt.errorbar(noms(x), noms(y), xerr=devs(x),yerr=devs(y),fmt='.', label='data')
plt.plot(smooth_x, fit.slope*smooth_x + fit.intercept, '--', label='fit')
plt.xlabel(r'x $\left[m\right]$')
plt.ylabel(r'$\frac{1}{\sqrt{R-R_b}}$ $\left[s^{1/2}\right]$')
plt.grid()
plt.legend()
plt.show()

# %%

thickness = np.array([0,    40,   80,   120,  160,  200, 240, 280, 320, 360, 400, 480, 560, ])*1e-4 # cm
time = np.array(     [20,   20,   30,   30,   40,   50,  80,  90,  100, 50,  50,  50 , 60,  ])
counts = np.array(   [1305, 1101, 1231, 1019, 1013, 1027,1233,1023,943, 389, 283, 151, 82])
# thickness = np.array([0,    40,   80,   120,  160,  200, 240, 280, 320, 360, 400, ])*1e-4 # cm
# time = np.array(     [20,   20,   30,   30,   40,   50,  80,  90,  100, 50,  50,    ])
# counts = np.array(   [1305, 1101, 1231, 1019, 1013, 1027,1233,1023,943, 389, 283,  ])

counts = uarray(counts, unp.sqrt(counts))
thickness = uarray(thickness, 1e-4)
rate = counts / time
normal_rate = rate - background_rate

density = 2.7e3 # mg/cm^3
density_thickness = thickness*density

plt.figure(dpi=130)
plt.errorbar(noms(density_thickness), noms(normal_rate), devs(normal_rate), devs(density_thickness), '.', label='data')
plt.xlabel(r'density thickness $\left[\frac{mg}{cm^2}\right]$')
plt.ylabel(r'cps $\left[\frac{1}{s}\right]$')
plt.grid()
plt.show()

#%%

eps = 0
d = (2*shelves_distance+top_shelf_GM_window+a)
y = unp.log(normal_rate[:-2]*d**2 + eps)
x = density_thickness[:-2]
fit = linregress(noms(x), noms(y))
mu = -ufloat(fit.slope, fit.stderr) # cm^2/mg
y = normal_rate[:-2]*d**2 + eps

plt.figure(dpi=130)
plt.errorbar(noms(x), noms(y), devs(y), devs(x), '.', label='data')
plt.plot(noms(x), np.exp(noms(x)*fit.slope+fit.intercept), '--', label='fit')
plt.xlabel(r'density thickness $\left[\frac{mg}{cm^2}\right]$')
plt.ylabel(r'cps $\left[\frac{1}{s}\right]$')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

range1 = -ufloat(fit.intercept,fit.intercept_stderr)/ufloat(fit.slope, fit.stderr)

E1 = unp.exp(6.63 - 3.2376*unp.sqrt(10.2146 - unp.log(range1)))

range2 = 8*np.log(2)/mu
E2 = unp.exp(6.63 - 3.2376*unp.sqrt(10.2146 - unp.log(range2)))


#%%

shelf_1 = 1141
shelf_1_time = 30
shelf_2 = 197
shelf_2_time = 30

disc = np.array([   0,  1,  2,  3,  4,  5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,])
time = np.array([   20, 20, 20, 20, 20, 20,20,20,20,20,20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,])
counts = np.array([ 10, 9,  11, 10, 11, 16,13,10,23,39,117,221,262,332,367,499,486,584,638,785,829])

dist0 = top_shelf_GM_window + 2*shelves_distance + a

disc_thick = 1e-3
x = dist0 - disc*disc_thick
counts = uarray(counts, unp.sqrt(counts))
rate = counts/time
normal_rate = (rate - background_rate) * (x+a)**2

curve_start = 9

fit = linregress(noms(x[curve_start:]), noms(normal_rate[curve_start:]))
smooth_x = np.linspace(np.min(noms(x)), -fit.intercept/fit.slope, 100)

mid = ufloat(noms(x)[10]+0.0003, 0.05e-2)
energy = ufloat(5.3,0.05) # MeV

plt.figure(dpi=130)
plt.errorbar(noms(x), noms(normal_rate), devs(normal_rate), devs(x), '.', label='data')
# plt.plot(smooth_x, smooth_x*fit.slope+fit.intercept, '--', label='fit')
plt.vlines(mid.n, 0, 0.08, 'darkorange', 'dashed', r'$R_{m}$')
# plt.ylim((-1,44))
plt.xlabel(r'd $\left[m\right]$')
plt.ylabel(r'cps $\cdot d^{2}$ $\left[\frac{m^2}{s}\right]$')
plt.legend()
plt.grid()
plt.show()


#%%
#Po-210
#sR90
# %%

a = np.linspace(0,6,100)
lmbda = 1
lmbda2 = 5
eps = -0.01
f = np.exp(-lmbda*a)
g = 0.5*(np.exp(-lmbda*a) + np.exp(-lmbda2*a))
plt.figure(dpi=130)
plt.plot(a, np.log(f), label='pure')
# plt.plot(a, np.log(np.exp(-lmbda2*a)), label='faster exp')
plt.plot(a, np.log(f+eps), label='neg const')
# plt.plot(a, np.log(f-eps), label='pos const')
# plt.plot(a, np.log(g), label='another exp')

plt.legend()
plt.grid()
plt.show()


