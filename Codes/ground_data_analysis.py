## Import libraries

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.stats

## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")

## Read input parameters

input_par_list = input_parameters.readlines()

bins = int(input_par_list[7])
y_max = float(input_par_list[8])

input_parameters.close()

## Define Guassina fit
## We set mu = 0 a priori

def Gaussian(x,Asquare,sigma):
    return Asquare*np.exp(-(x*x)/(sigma))

## Load data from simulations

y, dev_y = np.loadtxt(directory+'output/ground_bootstrap.txt', unpack=True)

interval = 2.0*y_max

x = np.linspace(-y_max + interval/(2*bins), y_max - interval/(2*bins), bins)

## Gaussian fit
init=(60.0,1.0)
popt, pcov = curve_fit(Gaussian, x, y, sigma=dev_y,p0=init)
chisq = (((y - Gaussian(x, *popt)) / dev_y)**2).sum()
ndof = len(y) - len(popt)
chisq = chisq/ndof

## Manipulate fit values and normalize fit function
sigma = popt[1]
dsigma = np.sqrt(pcov[1][1])

Integral = popt[0]*np.sqrt(popt[1]*np.pi)
dev_Integral = np.sqrt(pcov[0][0]*popt[1]*np.pi + popt[0]*popt[0]*pcov[1][1]*np.pi/(2*popt[1]))

A = popt[0]/Integral
dA = np.sqrt(pcov[0][0]/Integral**2 + popt[0]**2*dev_Integral/Integral**4)

A = np.sqrt(A)
dA = dA*0.5/A

print(f'A = {A} +- {dA}\nsigma = {sigma} +- {dsigma}\n')
print(f'chisq_norm = {chisq}')

## Plots

y = y/Integral # Plot normalized data points
dev_y = dev_y/Integral

plt.figure(1)

plt.grid(color = 'gray')
plt.title('Modulus squared of wave function')
plt.grid(color = 'gray',alpha=0.1)
plt.xlabel('position')
plt.errorbar(x,y,yerr=dev_y,color='plum', fmt='v')
plt.plot(x,Gaussian(x,A*A,sigma),color='blue')

plt.show()
