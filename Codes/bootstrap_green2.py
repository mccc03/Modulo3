## Import libraries

import string
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


## I/O directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")
output_green2=open(directory+'output/green2bootstrap.txt',"w")


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

y, y2 = np.loadtxt(directory+'means.txt', unpack=True, skiprows=thermal_len)
green2_data=np.loadtxt(directory+'green2.txt', unpack=True, skiprows=thermal_len)


mean = np.mean(y)
dev_mean = BootstrapDev(np.mean,y,resamplings)

for k in range(len(green2_data)):
    tmp_green = np.mean(green2_data[k]) - mean
    dev_tmp = BootstrapDev(np.mean,green2_data[k],resamplings)
    tmp_dev_green=np.sqrt(dev_tmp*dev_tmp+4*mean*mean*dev_mean*dev_mean)
    output_green2.write(str(tmp_green)+'\t'+str(tmp_dev_green)+'\n')

output_green2.close()
