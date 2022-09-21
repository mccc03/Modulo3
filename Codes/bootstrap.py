## Import libraries

import string
from tokenize import Double
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


## Initialize RNG

rng = np.random.default_rng()


## Define functions

def Energy(x,y):
    return 1.0/(2*eta)-x/(2*eta*eta)+y/2.0

def FitExp(x,A,l):
    return A*np.exp(-x/l)

def BootstrapDev(observable,array,resamplings):
    dev_list = []
    len_array = len(array)
    bin_size = int(len_array/16)
    while(bin_size <= 1+len_array/10):
        print('%d'%(bin_size))
        observable_bs = []
        for _ in range(resamplings):
            array_res = []
            for _ in range(int(len_array/bin_size)):
                i = rng.integers(len_array)
                array_res.extend(array[min(i,len_array-1-bin_size):min(i+bin_size,len_array-1)])
            observable_bs.append(observable(array_res))
        dev_list.append(np.std(observable_bs))
        bin_size*=2
    return max(dev_list)


## Load data from simulation

y, y2, delta_y2 = np.loadtxt(directory+'means.txt', unpack=True, skiprows=thermal_len)
green2_data=np.loadtxt(directory+'green2.txt', unpack=True, skiprows=thermal_len)
#green4=np.loadtxt(directory+'green4.txt', unpack=True, skiprows=thermal_len)

green2=[]
dev_green2=[]

mean = np.mean(y)
dev_mean = BootstrapDev(np.mean,y,resamplings)

for k in range(len(green2_data)):
    green2.append(abs(np.mean(green2_data[k]) - mean))
    dev_tmp = BootstrapDev(np.mean,green2_data[k],resamplings)
    dev_green2.append(np.sqrt(dev_tmp*dev_tmp+4*mean*mean*dev_mean*dev_mean))


#ene = Energy(mean_delta2,mean2)

x_fit = np.linspace(1.0, float(N/2), int(N/2))
x = np.linspace(1.0, float(N/2), 100)

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
