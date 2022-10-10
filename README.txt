INSTRUCTIONS TO RUN THE SIMULATIONS

sim_flag = 0

Run a simulation that calculates the energy of the system at a given temperature, the data is analyzed to verify the virial theorem and the behaviour of the kinetic energy. Steps:
 - compile HarmonicOscillator.cpp
 - run compiled cpp code
 - enter 0
 - run bootstrap_energy.py
 - repeat the previous steps for fixed T and different N's
 - run energy_data_analysis.py

Run a simulation to compute the energy at various temperatures and fit the data
 - compile HarmonicOscillator.cpp
 - run compiled cpp code
 - enter 0
 - run bootstrap_energy_temperature.py
 - repeat the previous steps for fixed eta and different T's
 - run energy_temperature_analysis.py

sim_flag = 1

Run a simulation that calculates the wave function of the ground state. Run this at low temperatures (T is actually the inverse temperature, so run at high T). Steps:
 - compile HarmonicOscillator.cpp
 - run compiled cpp code
 - enter 1
 - run bootstrap_ground.py
 - run ground_data_analysis.py

sim_flag = 2

Run a simulation that calculates the value of the first two energy gaps. Run this at low temperatures. Steps:
 - compile HarmonicOscillator.cpp
 - run compiled cpp code
 - enter 2
 - run bootstrap_green.py
 - repeat the previous steps for fixed T and different N's
 - run green2_data_analysis.py and green4_data_analysis.py
 - run gap1_fit.py and gap2_fit.py


The values in the file _data/input/parameters.txt correspond to the following variables:
 - T : inverse temperature of the system
 - measures : Number of measures taken
 - N : discretization depth
 - decorrel_len : distance between measures in units of Metropolis calls
 - thermal_len : guess for the number of steps required for the oscillator to reach equilibrium
 - init_flag : determines the starting state of the oscillator
 - resamplings : number of resamplings taken during bootstrap algorithm
 - bins: Number of bins for gaussian fit
 - y_max: for the gaussian fit, simulation bins lattice array from -y_max to y_max
