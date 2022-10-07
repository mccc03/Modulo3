## Import libraries

import string
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")
output_fit=open(directory+'output/green4_fit_T40.txt',"w")

output_fit.write('N\tA\tdA\t(E2-E0)\td(E2-E0)\tchi\n')

## Read input parameters

input_par_list = input_parameters.readlines()

T = float(input_par_list[0])

input_parameters.close()

## Define fit function

def FitExp(x,A,k):
    return A*np.exp(-x*k)


## Array of lattice dimensions

N_list = [400,286,222,182,154,133,118,105,95,87,80,74,69,65,61,57]

plt.figure(1)

plt.title('Fit for correlation length')
plt.ylabel(r'$\left\langle y_j y_{j+k}\right\rangle _c$')
plt.xlabel(r'$k$ [a]')
plt.grid(color = 'gray')

for N in N_list:
    ## Load data
    ## Data taken for T=1/40, N = 400/(1+0.4*n)

    green4,dev_green4=np.loadtxt(directory+'green4bootstrap_N'+str(N)+'.txt', unpack=True)

    ## Create linspace for fit and plot

    x_fit = np.linspace(1.0, float(len(green4)), int(len(green4)))
    x = np.linspace(1.0, float(len(green4)), int(len(green4)))

    ## Start fit

    popt, pcov = curve_fit(FitExp, x_fit, green4, sigma=dev_green4)
    chisq = (((green4 - FitExp(x_fit, *popt)) / dev_green4)**2).sum()
    ndof = len(green4) - len(popt)
    chisq = chisq/ndof

    ## Compute energy gap and propagate error (in units of \hbar\omega)
    eta = T/(float(N))

    gap = popt[1]/eta
    dev_gap = np.sqrt(pcov[1, 1])/eta

    ## Write fit results onto output file

    output_fit.write(str(N)+'\t'+str(popt[0])+'\t'+str(np.sqrt(pcov[0, 0]))+'\t'+str(gap)+'\t'+str(dev_gap)+'\t'+str(chisq)+'\n')

    ## Plot data points and best fit

    eta = T/float(N)

    plt.errorbar(x_fit,green4,yerr=dev_green4,fmt='.',label='eta = '+str(eta))
    plt.plot(x,FitExp(x,*popt),label='eta = '+str(eta))

output_fit.close()

plt.legend()
plt.show()
