Version 1.0

This set of python scripts runs our geological carbon cycle inverse model for the last 4.1 Ga. Posterior distributions for unknown parameters are plotted, along with time-evolution of carbon cycle variables. This version of the model is described in Krissansen-Totton, Kipp, and Catling (2021) "Carbon cycle inverse modeling suggests large changes in fractional organic burial are consistent with the carbon isotope record and may have contributed to the rise of oxygen", Geobiology. The carbon cycle model is adapted from Krissansen-Totton et al. (2018; PNAS), the original code for which is also available at github.com/joshuakt/early-earth-carbon-cycle.

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton et al. (2021). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (jkt@ucsc.edu).

REQUIREMENTS: Python, including numpy, pylab, scipy, emcee, and corner modules.

HOW TO RUN CODE:
(1) Put all the python scripts in the same directory, and ensure python is working in this directory.
(2) Open Main_code_inverse_model.py and check desired prior ranges, number of walkers, iterations, and threads (for parallelization)
(3) Run Main_code_inverse_model.py. When completed, outputs will be save and figures plotted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF CODE STRUCTURE:

%% Main_code_inverse_model.py
This script contains the likelihood function, and performs inverse calculations using emcee by calling the forward model. When inverse calculations are completed, it plots posterior distributions as corner plots, and calls another plotting script to show time-evolution of carbon cycle variables. Key input variables include:
- nwalk #Number of walkers in MCMC calculation
- ndim #Number of unknown paramteres
- nsteps  #Number of steps per walker, total forward model calls = nwalk*nsteps
- number_of_cores #Threads used to parallelize calculation 
Prior ranges for uncertain parameters can also be modified in the likelihood function, LnLike. 

%% forward_model
% Precambrian - Given parameter inputs, the forward model function Precambrian calculates the initial (Hadean) conditions for the carbon cycle e.g. equilibrium ocean chemistry and fluxes. Proportionality constants for carbon cycle functions are also calculated from initial conditions. Next, the ODE solver is called, and the system of equations describing the carbon cycle are solved. The ODE solver returns DIC and ALK as a function of time for both the ocean and the pore space. These outputs are fed back into carbon_cycle to obtain the time evolution for carbon cycle fluxes and ocean chemistry. Selected outputs are returned to Main_code_inverse_model.py. 

% system_of_equations - Contains the ODEs that describe the time evolution of the carbon cycle. The function takes the current state of the carbon and returns the time derivatives of DIC and ALK in the ocean and the pore space. This function is fed into the ODE solver to compute the time evolution of DIC, ALK for the ocean and pore space.

% lf - Calculates continental land fraction as a function of time.

% carbon_cycle - This function computes equilibrium chemistry and carbon cycle fluxes from DIC and ALK for both the ocean and the pore space. First, carbon_cycle calculates equilibrium chemistry for the ocean and the pore space, then heatflow and outgassing parameters. Next, it calculates surface and pore space temperatures from pCO2 using our climate model and surface-deep ocean relationship. Finally, carbon cycle fluxes are calculated from this information. Both carbon cycle fluxes and equilibrium chemistry are returned to system_of_equations or Precambrian.

%% plotting_everything.py
Contains plotting scripts for showing time-evolution of carbon cycle variables. 

## plotting_everything_for_paper.py
Similar to plotting_everything, except that after Main_code_invere_model has run once, this script can be used to produce figures without re-running forward model calculations. Additionally, this plotting script contains many of the figures from the main text in the Geobiology paper.

%% thermodynamic_variables.py
Contains the function Sol_prod that calculates carbonate solubility product as a function of temperature from Pilson, M. E. (1998) "An Introduction to the Chemistry of the Sea", Prentice-Hall, Inc. This function reproduces table G.1 in Pilson for salinity=35. thermodynamic_variables.py also contains the temperature-dependent functions for carbon equilibrium constants, equil_cont.

%% clim_model.py
This contains the the parameterized climate models used in this study.

%% Other files included in the code:
constraint_dict_revised.txt - contains true/observed values with uncertainties for comparison with model outputs in likelihood function
chain_output.npy - contains output chains from emcee calculations
lnprob_output.npy - contains corresponding log-likelihood values from emcee calculations
carb_200my_t.npy, carb_200my_t.npy, carb_200my_v.npy, carb_200my_var.npy - carbonate carbon isotope data
carb_RNG200my_t.npy, carb_RNG200my_t.npy, carb_RNG200my_v.npy, carb_RNG200my_var.npy - carbonate carbon isotope data
org_200my_t.npy,org_200my_t.npy, org_200my_v.npy, org_200my_var.npy - organic carbon isotope data
org_RNG200my_t.npy,org_RNG200my_t.npy, org_RNG200my_v.npy, org_RNG200my_var.npy - organic carbon isotope data
outputs_for_later_plotting3.npy - save time-evolution outputs for producing figures from paper using plotting_everything_for_paper script.

END EXPLANATION OF CODE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-------------------------
Contact e-mail: jkt@ucsc.edu
