#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import uarray
## Mo

#%% calibration part 1

# offset = 2%
# gain = level 2
# angle 1.9
# current 0.01 mA

#%% calibration part 2

def hist_peaks(file, height=40, prominence=30, verb=False):
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts #running average
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    Mes1 = pd.read_csv(file,sep='\t',header=1) # read the data.
    Counts = np.array(Mes1['Impulses/#']) # Impulses
    Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
    Channels = np.array(Mes1['Channel/#']) # Channel
    from scipy.signal import find_peaks
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
#%%
def fit_gauss1(channels, counts, p0):
    def Gauss(x, H, A, x0, sigma):
        return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    xdata=lin(channels,reg.slope,reg.intercept)
    ydata=counts

    parameters, covariance = curve_fit(Gauss, xdata, ydata,p0=p0)
    devs = np.sqrt(np.diag(covariance))

    H, A, x0, sigma = uarray(parameters, devs)
    
    fit_y = Gauss(xdata,  H.n, A.n, x0.n, sigma.n)

    plt.figure()
    plt.plot(xdata, ydata, 'o', label='data')
    plt.plot(xdata, fit_y, '-', label='fit')
    plt.grid()
    plt.legend()
    plt.show()

    return (H, A, x0, sigma)

#%% fit 2 gauss
def fit_gauss2(channels, counts, p0):
    def Gauss2(x, H, A1, x01, sigma1, A2, x02, sigma2):
        return H + A1 * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2)) + A2 * np.exp(-(x - x02) ** 2 / (2 * sigma2 ** 2))
    Channels=channels
    xdata=lin(Channels,reg.slope,reg.intercept)
    ydata=counts

    parameters, covariance = curve_fit(Gauss2, xdata, ydata, p0=p0)
    devs = np.sqrt(np.diag(covariance))
    
    H, A1, x01, sigma1, A2, x02, sigma2 = uarray(parameters, devs)
    
    fit_y = Gauss2(xdata, H.n, A1.n, x01.n, sigma1.n, A2.n, x02.n, sigma2.n)

    plt.figure(dpi=120)
    plt.plot(xdata, ydata, '.', label='data')
    plt.plot(xdata, fit_y, '-', label='fit')
    plt.xlabel('Channels')
    plt.ylabel('Impulses')
    plt.grid()
    plt.legend()
    return (H, A1, x01, sigma1, A2, x02, sigma2)

#%% Mo
peaks, chnls, counts, counts_smoothed = hist_peaks('data/Mo', verb=True)
plt.figure()
plt.plot(chnls[1900:2100], counts_smoothed[1900:2100])

#%% Zn
chnls = hist_peaks('data/Zn')

#%% Ni
chnls = hist_peaks('data/Ni')

#%% Ti
chnls = hist_peaks('data/Ti')

#%% Fe
chnls = hist_peaks('data/Fe', prominence=8)

#%% Cu
chnls = hist_peaks('data/Cu', height=15, prominence=10)

#%% Pb
chnls = hist_peaks('data/Pb', height=8, prominence=10)

#%%
# 42 Mo 17,479.34 17,374.3 19,608.3 
# 30 Zn 8,638.86 8,615.78 9,572.0
# 28 Ni 7,478.15 7,460.89 8,264.66
# 22 Ti 4,510.84 4,504.86 4,931.81
# 26 Fe 6,403.84 6,390.84 7,057.98
# 29 Cu 8,047.78 8,027.83 8,905.29
# 82 Pb 10,551.5 10,449.5 12,613.7 12,622.6 14,764.4
chan=np.array([2013,2267,937,1051,795,891,443,667,755,861,975,1170,1424,1696])
energies=np.array([17420, 19608,8625,9672,7470,8264,4648,6396,7058,8037,8905,10500,12618,14764])

reg=linregress(chan,energies)
xchan=np.arange(1,2500)
plt.figure(dpi=300)
plt.plot(chan ,energies,'.', label='Lines')
plt.plot(xchan, reg.slope*xchan + reg.intercept, '--', label='Fit')
plt.grid()
plt.ylabel('Energy [eV]')
plt.xlabel('Channels')
plt.legend()
plt.show()

a = ufloat(reg.slope, reg.stderr)
b = ufloat(reg.intercept, reg.intercept_stderr)

print(f'R squared: {reg.rvalue**2}')
#%% unknowns

def lin(x,a,b):
    return a*x+b
def e(x):
    return lin(x, a, b)

#%%unknown 1 (prob Cu)
#peaks=864, 968
#energies=8031.92388511, 8883.54218162
#8,047.78 8,027.83 8,905.29
chnls = hist_peaks('data/unknown1', prominence=10)
unknown1_E = e(chnls)

#%%unknown 2 (prob W)
#peaks=910, 1072, 1265
#e=8408.60120856,  9735.16047814, 11315.56750917
#8,397.6 8,335.2 9,672.35 9,961.5 11,285.9
chnls = hist_peaks('data/unknown2', height=3, prominence=5)
unknown2_E = e(chnls)


#%% crystal

# 27 Co 6,930.32 6,915.30 7,649.43
# 33 As 10,543.72 10,507.99 11,726.2

peaks, channels, counts, counts_smoothed = hist_peaks('data/CoAs', height=15, prominence=5, verb=True)
cryst_E = e(peaks)
plt.figure()
plt.plot(channels[1250:1375], counts_smoothed[1250:1375])
params1 = fit_gauss2(channels[600:800], counts_smoothed[600:800], (0,50,6400,200,25,7000,100))
params2 = fit_gauss1(channels[885:1000],counts_smoothed[885:1000], (0,20,8700,100))
# params3 = fit_gauss2
params3 = fit_gauss2(channels[1800:2100],counts_smoothed[1800:2100], (0,18,16900,300,20,17300,200))
# plt.figure()
# plt.plot(channels[1800:2100], counts_smoothed[1800:2100])

"""
6525.214591272148+/-44.497083359281454, - ???
6918.269189663606+/-45.23612252772048, - Co
8736.1467072241+/-49.1104433625483, - ???
10611.344687050014+/-53.7608782957335, - As
11700.433470093012+/-56.70379984498157 - As (low intensity)
16836.988412133465 dev 280 (2 gauss fit), prob Y
17403.076959542865 dev 210 (2 gauss fit), prob Mo
"""
"""
16,725.8 39 Y Kβ3 8
16,737.8 39 Y Kβ1 15
17,015.4 39 Y Kβ2 3
"""
"""
17,374.3 42 Mo Kα2 52
17,479.3 42 Mo Kα1 100
"""


"""
6,495.2 66 Dy 100
6,457.7 66 Dy 11
6,495.2 6,457.7 7,247.7 7,635.7 8,418.8
6,490.4 25 Mn 17
5,898.75 5,887.65 6,490.45
6,545.5 70 Yb 4
7,415.6 7,367.3 8,401.8 8,758.8 9,780.1
6,587.0 62 Sm 21
5,636.1 5,609.0 6,205.1 6,586 7,178
6,602.1 60 Nd 10
5,230.4 5,207.7 5,721.6 6,089.4 6,602.1
"""