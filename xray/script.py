import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
from scipy.optimize import curve_fit
## Mo

#%% calibration part 1

# offset = 2%
# gain = level 2
# angle 1.9
# current 0.01 mA

#%% calibration part 2

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts #running average
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
Mes1 = pd.read_csv('CoAs',sep='\t',header=1) # read the data.
Counts = np.array(Mes1['Impulses/#']) # Impulses
Counts_smoothed=smooth(Counts, 10) # smooth the data over 10 channels
Channels = np.array(Mes1['Channel/#']) # Channel
from scipy.signal import find_peaks
peaks, _ = find_peaks(Counts_smoothed, height=20, prominence=5)
plt.figure(dpi=300)
plt.plot(Channels ,Counts_smoothed, label='Original Mo-tube spectrum')
plt.plot(Channels[peaks] ,Counts_smoothed[peaks],"x",color='red', 
label='Lines')
plt.ylabel('Impulses')
plt.xlabel('Channels')
plt.legend()

#%%
# 42 Mo 17,479.34 17,374.3 19,608.3 
# 30 Zn 8,638.86 8,615.78 9,572.0
# 28 Ni 7,478.15 7,460.89 8,264.66
# 22 Ti 4,510.84 4,504.86 4,931.81
# 26 Fe 6,403.84 6,390.84 7,057.98
# 29 Cu 8,047.78 8,027.83 8,905.29
# 82 Pb 10,551.5 10,449.5 12,613.7 12,622.6 14,764.4
chan=np.array([2013,2267,937,1051,795,891,443,667,755,861,975,1170,1424,1696])
energies=np.array([17400, 19608,8625,9672,7470,8264,4648,6396,7058,8037,8905,10500,12618,14764])

reg=linregress(chan,energies)
xchan=np.arange(1,2500)
plt.figure(dpi=300)
plt.plot(chan ,energies,'.', label='Cal')
plt.plot(xchan, reg.slope*xchan + reg.intercept, '--')
plt.grid()
plt.ylabel('Energy [eV]')
plt.xlabel('Channels')
plt.legend()

#%% unknowns

def lin(x,a,b):
    return a*x+b
e=lin(peaks,reg.slope,reg.intercept)

#unknown 1 (prob Cu)
#peaks=864, 968
#energies=8031.92388511, 8883.54218162

#unknown 2 (prob W)
#peaks=910, 1072, 1265
#e=8408.60120856,  9735.16047814, 11315.56750917


#%% crystal

# 27 Co 6,930.32 6,915.30 7,649.43
# 33 As 10,543.72 10,507.99 11,726.2

#%% fit gauss
def Gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

xdata=lin(Channels[1250:1350],reg.slope,reg.intercept)
ydata=Counts_smoothed[1250:1350]

parameters, covariance = curve_fit(Gauss, xdata, ydata,p0=(5,120,10600,500))
  
H, A, x0, sigma = parameters
  
fit_y = Gauss(xdata,  H, A, x0, sigma)

plt.plot(xdata, ydata, 'o', label='data')
plt.plot(xdata, fit_y, '-', label='fit')
plt.legend()

#%% fit 2 gauss

def Gauss2(x, H, A1, x01, sigma1, A2, x02, sigma2):
    return H + A1 * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2)) + A2 * np.exp(-(x - x02) ** 2 / (2 * sigma2 ** 2))

xdata=lin(Channels[1250:1350],reg.slope,reg.intercept)
ydata=Counts_smoothed[1250:1350]

parameters, covariance = curve_fit(Gauss2, xdata, ydata,p0=(5,25,11700,300,5,11950,100))
  
H, A1, x01, sigma1, A2, x02, sigma2 = parameters
  
fit_y = Gauss2(xdata, H, A1, x01, sigma1, A2, x02, sigma2)

plt.plot(xdata, ydata, 'o', label='data')
plt.plot(xdata, fit_y, '-', label='fit')
plt.legend()