## Import libraries

import string
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl

## Random color in plot

import random


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")
output_fit=open(directory+'output/green2_fit_T40.txt',"w")

output_fit.write('N\tA\tdA\t(E1-E0)\td(E1-E0)\tchi\n')

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
plt.xlabel('$k$ [a]')
plt.grid(color = 'gray',alpha=0.4)

#colors = plt.cm.jet(np.linspace(0, 1, len(N_list)))
cmap=mpl.colormaps['turbo']
colors = cmap(N_list)

i = 0

for N in N_list:
    ## Load data
    ## Data taken for T=1/40, N = 400/(1+0.4*n)

    green2,dev_green2=np.loadtxt(directory+'green2bootstrap_N'+str(N)+'.txt', unpack=True)

    ## Create linspace for fit and plot

    x_fit = np.linspace(1.0, float(len(green2)), int(len(green2)))
    x = np.linspace(1.0, float(len(green2)), int(len(green2)))

    ## Start fit

    popt, pcov = curve_fit(FitExp, x_fit, green2, sigma=dev_green2)
    chisq = (((green2 - FitExp(x_fit, *popt)) / dev_green2)**2).sum()
    ndof = len(green2) - len(popt)
    chisq = chisq/ndof

    ## Compute energy gap and propagate error
    eta = T/(float(N))

    gap = popt[1]/eta
    dev_gap = np.sqrt(pcov[1, 1])/eta

    ## Write fit results onto output file

    output_fit.write(str(N)+'\t'+str(popt[0])+'\t'+str(np.sqrt(pcov[0, 0]))+'\t'+str(gap)+'\t'+str(dev_gap)+'\t'+str(chisq)+'\n')

    ## Plot data points and best fit

    eta = T/float(N)

    color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

    plt.errorbar(x_fit,green2,yerr=dev_green2,fmt='.',label=rf'$\eta =${eta:.3f}',color=colors[i])
    plt.plot(x,FitExp(x,*popt),color=colors[i])

    i = i+1

output_fit.close()

plt.legend()
plt.show()
