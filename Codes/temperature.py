## Import libraries

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## Define functions

def Energy(x,y):
    return 5.0 - 50.0*x + 0.5*y

def FitEne(x,a,b):
    return a - b/(np.exp(1/x)-1)


## Set working directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/output/'


## Read data

T,y,dy,y2,dy2,Dy2,dDy2 = np.loadtxt(directory+'energy_temp.txt', unpack=True)

## Compute energy and error

ene = Energy(Dy2,y2)
dev_ene = 50.0*dDy2

## Start fit
initial=(0.5,1.0)
popt, pcov = curve_fit(FitEne, T, ene, sigma=dev_ene,absolute_sigma=False,p0=initial)
chisq = (((ene - FitEne(T, *popt)) / dev_ene)**2).sum()

ndof = len(T)-len(popt)

print(f'a = {popt[0]:.4f} +/- {np.sqrt(pcov[0, 0]):.4f}')
print(f'b = {popt[1]:.4f} +/- {np.sqrt(pcov[1, 1]):.4f}')
print(f'norm chisq = {chisq/ndof}')

plt.figure(1)
plt.title('Energy as a function of temperature')
plt.xlabel(r'$T = 1/(N\eta)$')
plt.ylabel('Energy')
plt.grid(color = 'gray')
plt.errorbar(T,ene,yerr=dev_ene, color='green',fmt='.')
plt.plot(T,FitEne(T,*popt), color='blue')

plt.show()


