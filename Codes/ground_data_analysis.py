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

measures = int(input_par_list[0])
resamplings = int(input_par_list[7])
eta = float(input_par_list[2])
N = int(input_par_list[1])
thermal_len = int(measures*N/2)

input_parameters.close()

## Define function

def Gaussian(x,A,mu,sigma):
    return A*np.exp(-(x-mu)*(x-mu)/(sigma))

## Load data from simulations

y = np.loadtxt(directory+'position.txt', unpack=True, skiprows=thermal_len)


## Gaussian fit

bins_number = 50
y_norm, bins, _ = plt.hist(y,bins=bins_number,density=True)
bins1=np.zeros(bins_number)
for i in range(bins_number):
    bins1[i] = 0.5*(bins[i+1]+bins[i])


x_fit = np.linspace(y.min(), y.max(), 500)
popt, pcov = curve_fit(Gaussian, bins1, y_norm)
print(f'A = {popt[0]} +- {np.sqrt(pcov[0][0])}\nmu = {popt[1]} +- {np.sqrt(pcov[1][1])}\nsigma = {popt[2]} +- {np.sqrt(pcov[2][2])}\n')

## Plots

plt.figure(1)

plt.grid(color = 'gray')
plt.title('Modulus squared of wave function')
plt.grid(color = 'gray')
plt.xlabel('position')
plt.hist(y,bins=bins_number,density=True,color='green')
plt.plot(x_fit,Gaussian(x_fit,*popt),color='blue')

plt.show()
