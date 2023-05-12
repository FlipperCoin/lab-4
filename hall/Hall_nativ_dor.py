#%% init

import dis
import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from scipy.optimize import curve_fit as cfit
from scipy.constants import k as Kb
from scipy.constants import eV

from uncertainties.unumpy import uarray
from uncertainties import unumpy
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import nominal_values as uval
from uncertainties.unumpy import std_devs as uerr


# helper function for plotting data and regression
def one4all(xdata,ydata,yerr=0,xerr=0,mode="none",f=None,xlabel="x",ylabel="y",show=True,label="Data"):
    fig = plt.figure(dpi=120)
    plt.errorbar(xdata,ydata,yerr,xerr,".",label=label)

    
    if mode == "none":
        fit= []
        
        
    
    elif mode == "linear":
        fit = linregress(xdata,ydata)
        f = lambda x,a,b: a*x+b
        #fit_label = "Regression: y=" + str(fit.slope) + str("x+") + str(fit.intercept)
        plt.plot(xdata,f(xdata,fit.slope,fit.intercept),"-.",label="Regression")
    
       
        
    elif mode == "0 intercept":
        f = lambda x,a: a*x
        fit = cfit(f,xdata,ydata)
        plt.plot(xdata,f(xdata,*fit[0]),"-.",label="Regression")
        
    elif mode == "general function":
        if f==None:
            raise TypeError

        fit = cfit(f,xdata,ydata)
        plt.plot(xdata,f(xdata,*fit[0]),"-.",label="fit")
    
    
    else:
        print("mode='",mode,"' is not definted!")
        raise TypeError


    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.grid()
    if mode != "none":
        plt.legend()
        
    if show:
        plt.show()
    
    return (fig,fit)

def Rsqrue(x,y):
    RSS = np.sum((y-x)**2)
    TSS = np.sum(y**2)
    return 1 - RSS/TSS

def Reg_print(fit):
    m = ufloat(fit.slope,fit.stderr)
    b = ufloat(fit.intercept,fit.intercept_stderr)
    print("==> y =(",m,")x + (",b,") , R^2=",fit.rvalue**2)

#%% part 0
                            ######### Station Number: 3 #########

d = 1e-3 #m
L = 16e-3 #m
W = 10e-3 #m
e = 1.6e-19 #C

Ip = np.array([    -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30]) * 1e-3 #A
I_err = 1e-3 #A
Up = np.array([     -1.668, -1.420,  -1.127,  -0.901,  -0.628, -0.393,   -0.133, 0.217, 0.433,    0.739,  1.019,  1.233,  1.495]) #V
U_err = 0.001 #V

fig,fit = one4all(xdata=Ip,ydata=Up,xerr=I_err,yerr=U_err,xlabel=r"$I_p$ $\left[A\right]$",ylabel=r"$U_p$ $\left[V\right]$",mode="linear")
Reg_print(fit)

R0 = ufloat(fit.slope,fit.stderr)
print("r=",R0,"ohm")

rho0 = R0*d*W/L

#%% part 1

B = 251e-3 #T
B_err =1e-3 #T
B = ufloat(B,B_err)

Ip = np.array([     -30,     -25,   -20,     -15,    -10,   -5,     0,   5,   10,     15,     20,     25,     30]) * 1e-3 #A
I_err = 1e-3 #A
UH = np.array([     60.00,   49.39,  39.88,   31.70,  22.99,12.92,1.99,-8.86,-16.44, -26.26, -35.87, -45.19,   -54.67])*1e-3 #V
U_err = 1e-4 #V

fig,fit = one4all(xdata=Ip,ydata=UH,xerr=I_err,yerr=U_err,xlabel=r"$I_p[A]$",ylabel=r"$U_H[V]$",mode="linear")
Reg_print(fit)

slope = ufloat(fit.slope,fit.stderr)
RH1 = slope * d / B
print("RH=",RH1,"m^-3*C-1")

if RH1<0:
    print("The s/c is of P-type!")
else:
    print("The s/c is of N-type!")

n1 = 1/np.abs(RH1)/e
print("n=",n1,"m^-3")

mu1 = np.abs(RH1)/rho0
print("mu=",mu1,"m^-4*(ohm*C)-1")


#%% part 2
I = 30e-3 #A
I_err = 1e-3
I = ufloat(I,I_err)

B = np.array([    -300, -280, -260, -240, -220, -200, -180, -160, -140, -120, 
                  -100,  -80,  -60,  -40,  -20,    0,   20,   40,   60,   80,
                   100,  120,  140,  160,  180,  200,  220,  240,  260,  280,
                   300]) * 1e-3 #T
B_err = 1e-3 #T
#we took the measures from 300 to -300
UH = np.array([    -61.52,  -57.90, -54.29, -50.38, -46.57, -42.50, -38.59, -34.42, -30.13,  -25.95,    -21.32, -17.13, -12.65, -7.84,  -3.15,  1.60,
                   6.16,    10.87,  15.37,  19.97,  24.38,  28.79,  33.06,  37.38,  41.59,  45.56,      49.49,  53.44,  57.15,  60.93,  64.64])[::-1]*1e-3 #V
U_err = 1e-5 #V


fig,fit = one4all(xdata=B,ydata=UH,xerr=B_err,yerr=U_err,xlabel=r"$B[T]$",ylabel=r"$U_H[V]$",mode="linear")
Reg_print(fit)

slope = ufloat(fit.slope,fit.stderr)
RH2 = slope * d / I
print("RH=",RH2,"m^-3*C-1")

if RH2<0:
    print("The s/c is of P-type!")
else:
    print("The s/c is of N-type!")

n2 = 1/np.abs(RH2)/e
print("n=",n2,"m^-3")

mu2 = np.abs(RH2)/rho0
print("mu=",mu2,"m^-4*(ohm*C)-1")

rel_dif = np.abs(RH1-RH2)/max(np.abs(RH1.n),np.abs(RH2.n)) * 100
print("The differnce between RH1 to RH2 is:",rel_dif,"%")

#%% part 3
Ip = 30e-3 # A
Ip_err = 0.1e-3 #A
Ip = ufloat(Ip,Ip_err)

# Take only the higest precesion
# There was an offset of 7 mT on the first measurements duo to residule magnetic field
B = np.array([     6    ,10      ,30       ,50     ,70      ,90     ,100       ,120     ,140        ,160        ,180        ,200        ,220        ,240        ,260        ,280        ,300]) * 1e-3 #T
Up = np.array([    1.560,1.5602  ,1.5608   ,1.562  ,1.5632  ,1.5647 ,1.5660    ,1.5683  ,1.571      ,1.574      ,1.5773     ,1.5808     ,1.5846     ,1.5887     ,1.5927     ,1.5971     ,1.6017]) #V
B_err = 1e-3#T
Up_err = 1e-4 #V

B = uarray(B,B_err)
Up = uarray(Up,Up_err)

R = Up/Ip

x = B**2
y = R - R0
fig,fit = one4all(xdata=uval(x),ydata=uval(y),xerr=uerr(x),yerr=uerr(y),xlabel=r"$B^{2}\, \left[T^{2}\right]$",ylabel=r"$\left(R - R_{0}\right)\,[\Omega]$",mode="linear",show=True)
# axes = fig.axes
# axes[0].set_ylim(30, 60)
# fig.show()
Reg_print(fit)

def f(x, a, b):
    return a + b*x**2

fig,fit = one4all(xdata=uval(B),ydata=uval(R),xerr=uerr(B),yerr=uerr(R),xlabel=r"$B\, \left[T\right]$",ylabel=r"$R\,[\Omega]$",mode="general function",f=f, show=True)
R0_2 = ufloat(fit[0][0], np.sqrt(fit[1][0][0]))

mu_e_ref = ufloat(0.39, 0.01)
ni = 2.4*10**13 / (0.1)**3
p = n1
n_e = ni**2/p
rho0_theo = 1/(e*(mu_e_ref*n_e+mu1*p))
R0_theo = rho0_theo*L/W/d
rho_xx = R0_2*W*d/L + (B**2) * mu_e_ref * n_e * mu1 * p * (mu_e_ref+mu1)**2 / (e * (mu_e_ref*n_e + mu1*p)**3)
R_2carrier = rho_xx*L/W/d
rho_xy = B*(mu1**2*p-mu_e_ref**2*n_e)/(e*(mu_e_ref*n_e+mu1*p)**2)
theta = rho_xy/rho_xx
R_geo = R0_2 + 4*rho_xx/(np.pi*d) * theta**2 / ((np.pi/2)**2 - theta**2)

fig,fit = one4all(xdata=uval(B),ydata=uval(R),xerr=uerr(B),yerr=uerr(R),xlabel=r"$B\, \left[T\right]$",ylabel=r"$R\,[\Omega]$",mode="general function",f=f, show=False)
plt.figure(fig)
plt.errorbar(uval(B), uval(R_2carrier), xerr=uerr(B), yerr=uerr(R_2carrier),label='2 carriers')
plt.errorbar(uval(B), uval(R_geo), xerr=uerr(B), yerr=uerr(R_geo),label='geometric')
plt.legend()
plt.show


#%% part 4
I = 30e-3 #A
I_err = 1e-3
I = ufloat(I,I_err)

T = np.append(np.arange(30,111, 10), np.arange(120,141,1)) + 273.15 #K
T_err = np.array([ 0]) #C
Up = np.append(np.array([1.5889, 1.6950, 1.7942, 1.8710, 1.8837, 1.7948, 1.5679, 1.2828, 0.9825]), np.array([    7343, 7068, 6876, 6718, 6549, 6299, 6156, 5884, 5771, 5659, 5455, 5319, 5201, 5073, 4895, 4727, 4637, 4599, 4470, 4326, 4211])*1e-4) #V
Up_err = np.array([ 0]) #V

T = uarray(T, 1)
Up = uarray(Up, 0.0001)
# T = uarray(T,T_err)
# Up = uarray(Up,Up_err)

fig,fit = one4all(xdata=uval(1/T),ydata=uval(1/Up),xerr=uerr(1/T),yerr=uerr(1/Up), xlabel=r"$T^{-1}\,\left[K^{-1}\right]$",ylabel=r"$U_p^{-1}\,\left[V^{-1}\right]$",mode="none")

fig,fit = one4all(xdata=uval(1/T),ydata=uval(-unumpy.log(Up)),xerr=uerr(1/T),yerr=uerr(-unumpy.log(Up)),xlabel=r"$T^{-1}\,\left[K^{-1}\right]$",ylabel=r"-$\ln\left(U_p\right)$",mode="none",show=False)
x = 1/T
y = -unumpy.log(Up)
y = y[x < 0.0026]
x = x[x < 0.0026]
fit = linregress(uval(x), uval(y))
x_smooth_max = (np.min(-unumpy.log(Up)) - fit.intercept) / fit.slope
x_smooth = np.linspace(np.min(x).n, x_smooth_max.n,100)

plt.figure(fig)
plt.plot(x_smooth, fit.intercept + fit.slope*x_smooth, '--', label='fit')
plt.legend()
plt.show()

Eg1 = ufloat(fit.slope, fit.stderr)*(-2*Kb) / eV, 

#%% part 5
B= 300e-3 #T
B_err =1e-3 
B = ufloat(B,B_err)

I = 30e-3 #A
I_err = 1e-3
I = ufloat(I,I_err)

TUH = [(137,4.94), (136,5.02),(134,5.04),  (133,5.14), (132,5.20), (131,5.25), (130,5.30), (129,5.32), (128,5.37), (127,5.41), (126,5.45), (125,5.47),(123,5.47), (120,5.41),
       (119,5.36),  (117,5.24), (115,4.95), (112,4.4),(109,3.24),(106,1.99),    (102,-0.59),(99,-3.48),(93,-11),(85,-23.23),    (75,-40.71),(65,-53.75),
       (55,-60.52), (45,-62.96),    (35,-63.43),    (31,-63.27)]
T = []
UH =[]
for t,uh in TUH:
    T.append(t)
    UH.append(uh)

T = np.array(T) + 273.15 #K
UH = np.array(UH)*1e-3 #V

# UH = UH[T>(130 + 273.15)]
# T = T[T>(130 + 273.15)]

T_err = 1 #k
UH_err = 1e-4 #V

T = uarray(T,T_err)
UH = uarray(np.abs(UH),UH_err)

x = 1/T
y = unumpy.log(UH)
fig,fit = one4all(xdata=uval(x),ydata=uval(y),xerr=uerr(x),yerr=uerr(y),xlabel=r"$T^{-1}\left[K^{-1}\right]$",ylabel=r"$\ln\left(\left|U_{H}\right| \right)$",mode="none",show=False)
y = y[:12]
x = x[:12]
smooth_x = np.linspace(np.min(x).n, 0.0026, 100)
fit = linregress(uval(x), uval(y))
plt.figure(fig)
plt.plot(smooth_x, fit.slope*smooth_x + fit.intercept, '--', label='fit')
plt.legend()
plt.show()

Eg2 = ufloat(fit.slope, fit.stderr)*(2*Kb) / eV, 