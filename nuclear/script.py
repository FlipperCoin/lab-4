# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# %% [markdown]
# ### Statistics of counting and Background Radiation Measurement

# %%
background_count = 0
background_rate = background_count / 100 # per sec

# %%
init_nbar = np.sum([ ]) / 10
m = 150*init_nbar

# %%
df = pd.read_csv('data/data.csv')

# %%
n = df['n']

# %%
print(f"m = {m} and len(n) = {len(n)}")

nbar = np.mean(n)
n_std = np.sqrt(m) * np.sqrt(np.mean((n-nbar)**2))

K3 = (n-nbar)**3
K3_bar = np.mean(K3)
K3_std = np.sqrt(m) * np.sqrt(np.mean((K3-K3_bar)**2))

# %%
print(f"n +/- STD(nbar) = {nbar} +/- {n_std}")
print(f"K3 +/- STD(K3) = {K3_bar} +/- {K3_std}")

# %% [markdown]
# ### Inverse Square Law

# %%
shelves_distance = 10e-3
top_shelf_GM_window = 12.3e-3
stage_thickness = 1e-3
number_of_shelves = 5

time = 0
count_per_shelf = np.array([1, 2, 3, 4, 5, 6])

rate_per_shelf = count_per_shelf / time
normal_rate_per_shelf = rate_per_shelf - background_rate

distance_per_shelf = number_of_shelves*(shelves_distance-1) + top_shelf_GM_window - np.arange(number_of_shelves)*shelves_distance

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


