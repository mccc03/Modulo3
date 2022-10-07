## Import libraries

import string
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")

## Read input parameters

input_par_list = input_parameters.readlines()

T = float(input_par_list[0])

input_parameters.close()

## Define fit function

def LinearFit(x,m,q):
    return m*x + q

## Array of lattice dimensions

N_list = [400,286,222,182,154,133,118,105,95,87,80,74,69,65,61,57]

## Create fit data

## X axis --> eta^2

eta2_list = []
for N in N_list:
    eta2_list.append(T*T/(float(N)*float(N)))
eta2=np.array(eta2_list)

## Y axis --> gap

N,A,dev_A,gap,dev_gap,chi = np.loadtxt(directory+'output/green4_fit_T40.txt', unpack=True, skiprows=1)

## Start fit

popt, pcov = curve_fit(LinearFit, eta2, gap, sigma=dev_gap)
chisq = (((gap - LinearFit(eta2, *popt)) / dev_gap)**2).sum()
ndof = len(gap) - len(popt)
chisq = chisq/ndof

## Print fit results

print(f'm = {popt[0]} +- {np.sqrt(pcov[0][0])}\nE1-E0 = {popt[1]} +- {np.sqrt(pcov[1][1])}')
print(f'chisq_norm = {chisq}')


## Plot results

plt.figure(1)

plt.title('Fit for second energy gap')
plt.ylabel(r'$E2-E0 [\hbar\omega]$')
plt.xlabel(r'$\eta^2$')
plt.grid(color = 'gray')
plt.errorbar(eta2,gap,yerr=dev_gap,fmt='.')
plt.plot(eta2,LinearFit(eta2,*popt))

plt.show()
