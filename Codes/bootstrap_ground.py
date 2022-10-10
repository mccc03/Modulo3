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

input_parameters.close()

## Set output files

output_ground=open(directory+'output/ground_bootstrap.txt', "w")

## Initialize RNG

rng = np.random.default_rng()

## Define functions

def BootstrapDev(observable,array,resamplings):
    dev_list = []
    len_array = len(array)
    bin_size = int(len_array/(1024))
    while(bin_size <= 1+len_array/10):
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

prob_data = np.loadtxt(directory+'probs.txt', unpack=True, skiprows=thermal_len)

## Call bootstrap for each bin

for k in range(len(prob_data)):
    tmp_prob = np.mean(prob_data[k])
    tmp_dev = BootstrapDev(np.mean,prob_data[k],resamplings)
    output_ground.write(str(tmp_prob)+'\t'+str(tmp_dev)+'\n')

