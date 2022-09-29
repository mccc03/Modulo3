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

thermal_len = int(input_par_list[4])
resamplings = int(input_par_list[7])
eta = float(input_par_list[2])
N = int(input_par_list[1])

input_parameters.close()

## Define functions

def FitExp(x,A,l):
    return A*np.exp(-x/l)


green2,dev_green2=np.loadtxt(directory+'output/green2bootstrap.txt', unpack=True)

x_fit = np.linspace(1.0, float(len(green2)), int(len(green2)))
x = np.linspace(1.0, float(len(green2)), 100)

popt, pcov = curve_fit(FitExp, x_fit, green2, sigma=dev_green2)
chisq = (((green2 - FitExp(x_fit, *popt)) / dev_green2)**2).sum()

print(f'A = {popt[0]:.4f} +/- {np.sqrt(pcov[0, 0]):.4f}')
print(f'l = {popt[1]:.4f} +/- {np.sqrt(pcov[1, 1]):.4f}')
print(f'cov = {np.sqrt(pcov[1, 0]):.4f}')

plt.figure(1)

plt.title('Fit for first energy gap')
plt.ylabel('2 point connected Green function')
plt.xlabel('Distance')
plt.grid(color = 'gray')
plt.errorbar(x_fit,green2,yerr=dev_green2, color='green',fmt='.')
plt.plot(x,FitExp(x,popt[0],popt[1]), color='blue')

plt.show()
