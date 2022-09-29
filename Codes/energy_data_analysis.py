## Import libraries

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## Define functions

def Energy(x,y,k):
    return 1.0/(2*k) - x/(2*k*k) + y/2.0

def KineticNorm(x,k):
    return 1.0/(2*k) - x/(2*k*k)

def KineticDiv(x,k):
    return - x/(2*k*k)   

def Parabola(x,a,b):
    return a*x*x + b

## Set working directory

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/output/'

## Read data

N,y,dy,y2,dy2,Dy2,dDy2 = np.loadtxt(directory+'energy20.txt', unpack=True, skiprows=2)

## Create eta variable

eta=20.0/N

## Start fit

popt, pcov = curve_fit(Parabola, eta, y2, sigma=dy2,absolute_sigma=False)
chisq = (((y2 - Parabola(eta, *popt)) / dy2)**2).sum()

ndof = len(y2)-len(popt)

print(f'a = {popt[0]:.4f} +/- {np.sqrt(pcov[0, 0]):.4f}')
print(f'b = {popt[1]:.4f} +/- {np.sqrt(pcov[1, 1]):.4f}')
print(f'norm chisq = {chisq/ndof}')

x=np.linspace(0.0,0.5,100)

plt.figure(1)

plt.title('Potential energy')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\left\langle y^2\right\rangle$')
plt.grid(color = 'gray')
plt.errorbar(eta,y2,yerr=dy2, color='green',fmt='.')
plt.plot(x,Parabola(x,*popt), color='blue')

plt.figure(2)

plt.title('Divergent kinetic energy')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$-\frac{\left\langle \Delta y^2\right\rangle}{2\eta ^2}$')
plt.grid(color = 'gray')
plt.errorbar(eta,KineticDiv(Dy2,eta),yerr=(dDy2/(2*eta*eta)), color='green',fmt='.')


plt.figure(3)
plt.title('Normalized kinetic energy')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\frac{1}{2\eta} - \frac{\left\langle \Delta y^2\right\rangle}{2\eta ^2}$')
plt.grid(color = 'gray')
plt.errorbar(eta,KineticNorm(Dy2,eta),yerr=(dDy2/(2*eta*eta)), color='blue',fmt='.')

plt.show()
