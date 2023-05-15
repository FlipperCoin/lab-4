# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# %% [markdown]
# ### Statistics of counting and Background Radiation Measurement

# %%
# range 500 - 1200 v, steps 15v, 3 seconds per step, 
plat_df = pd.read_csv('D:/itai/data/plat.tsv', sep='\t+', header=9)

plt.plot(plat_df['Voltage'], plat_df['Counts'], '.')
plt.grid()

op_v = 950 # v

# %%
background_count = 43
background_rate = background_count / 100 # per sec

# %%
init_nbar = np.sum([0,5,5,5,2,3,3,1,4,2]) / 10
m = 150*init_nbar

# %%
df = pd.read_csv('D:/itai/data/nbar3.tsv', sep='\t+', header=9)

# %%
n = df['Counts']

# %%
print(f"m = {m} and len(n) = {len(n)}")

nbar = np.mean(n)
n_std = np.sqrt(np.mean((n-nbar)**2)) # *np.sqrt(m)

K3 = (n-nbar)**3
K3_bar = np.mean(K3)
K3_std = np.sqrt(np.mean((K3-K3_bar)**2)) # *np.sqrt(m)

print(f"n +/- STD(nbar) = {nbar} +/- {n_std}")
print(f"K3 +/- STD(K3) = {K3_bar} +/- {K3_std}")


# %% [markdown]
# ### Inverse Square Law

# %%
shelves_distance = 1.09e-2 # guide ~10e-3
top_shelf_GM_window = 1.09e-2 # guide 12.3e-3
# from second shelf 1.94e-3 
stage_thickness = 1e-3
number_of_shelves = 10

time_per_shelf = np.array([ 60,  50,    50,   40,    30,    20,    20,  8,    5,   3])
count_per_shelf = np.array([1320,1198 ,1546, 1527, 1493,  1312,   1936,  1154, 1209,1303])

rate_per_shelf = count_per_shelf / time_per_shelf
normal_rate_per_shelf = rate_per_shelf - background_rate

distance_per_shelf = number_of_shelves*(shelves_distance-1) + top_shelf_GM_window - np.arange(number_of_shelves)*shelves_distance

plt.plot(rate_per_shelf, distance_per_shelf, '.')
plt.show()

x = distance_per_shelf
y = 1/np.sqrt(rate_per_shelf)
fit = linregress(x, y)

a = fit.intercept/fit.slope

# %%
smooth_x = np.linspace(np.min(x), np.max(x), 100)
plt.plot(x, y, '.', label='data')
plt.plot(smooth_x, fit.slope*smooth_x + fit.intercept, '--', label='fit')
plt.grid()
plt.legend()
plt.show()

# %%

thickness = np.array([0,    40,   80,   120,  160,  200, 240, 280, 320, 360, 400, 480, 560, ])*1e-4 # cm
time = np.array(     [20,   20,   30,   30,   40,   50,  80,  90,  100, 50,  50,  50 , 60,  ])
counts = np.array(   [1305, 1101, 1231, 1019, 1013, 1027,1233,1023,943, 389, 283, 151, 82])

rate = counts / time
normal_rate = rate - background_rate

density = 2.7e3 # mg/cm^3
density_thickness = thickness*density

plt.plot(density_thickness, normal_rate, '.', label='data')
plt.grid()
plt.show()

#%%

y = np.log10(normal_rate)
x = density_thickness
fit = linregress(x, y)
mu = -fit.slope

plt.plot(x, y, '.', label='data')
plt.plot(x, x*fit.slope+fit.intercept, '--', label='fit')
plt.grid()
plt.show()


#%%

shelf_1 = 1141
shelf_1_time = 30
shelf_2 = 197
shelf_2_time = 30

disc = np.array([   0,  1,  2,  3,  4,  5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,])
time = np.array([   20, 20, 20, 20, 20, 20,20,20,20,20,20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,])
counts = np.array([ 10, 9,  11, 10, 11, 16,13,10,23,39,117,221,262,332,367,499,486,584,638,785,829])

dist0 = top_shelf_GM_window + 2*shelves_distance

disc_thick = 1e-3
x = dist0 - disc*disc_thick
rate = counts/time
normal_rate = rate - background_rate

plt.plot(x, normal_rate, '.')
plt.grid()
plt.show()

#%%
#Po-210
#sR90