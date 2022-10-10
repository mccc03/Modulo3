## Import libraries

import string
import numpy as np


## Set input directory and files

directory = '/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/'
input_parameters=open(directory+'input/parameters.txt',"r")


## Read input parameters

input_par_list = input_parameters.readlines()

thermal_len = int(input_par_list[4])
resamplings = int(input_par_list[6])
N = int(input_par_list[2])

input_parameters.close()


## Set output files

output_energy=open(directory+'output/energy_temp.txt', "a")


## Initialize RNG

rng = np.random.default_rng()


## Define functions

def BootstrapDev(observable,array,resamplings):
    dev_list = []
    len_array = len(array)
    bin_size = int(len_array/128)
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

mean = np.mean(y)
dev_mean = BootstrapDev(np.mean,y,resamplings)
mean2 = np.mean(y2)
dev_mean2 = BootstrapDev(np.mean,y2,resamplings)
mean_delta2 = np.mean(delta_y2)
dev_mean_delta2 = BootstrapDev(np.mean,delta_y2,resamplings)

## temperature = 1.0/(eta*N) = 10.0/N
## eta = 0.1

temperature = float(10.0/float(N))

output_energy.write(str(temperature)+'\t'+str(mean)+'\t'+str(dev_mean)+'\t'+str(mean2)+'\t'+str(dev_mean2)+'\t'+str(mean_delta2)+'\t'+str(dev_mean_delta2)+'\n')

output_energy.close()
