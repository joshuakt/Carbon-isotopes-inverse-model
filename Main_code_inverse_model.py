#########################################
## Load required modules
import emcee
import numpy
import matplotlib.pyplot as pl
from corner import corner
import corner
from forward_model import Precambrian ## The function 'Precambrian' is the forward model
from scipy.interpolate import interp1d
#########################################


#########################################
## Load binned carbonate carbon isotope data
carb_filt=numpy.ones(shape=(len(numpy.load("carb_200my_t.npy")),3))
carb_filt[:,0]=numpy.load("carb_200my_t.npy") #time values for carbonates
carb_filt[:,1]=numpy.load("carb_200my_v.npy") #isotope values for carbonates
carb_filt[:,2]=numpy.load("carb_200my_var.npy") #isotope values for carbonates

carb_filt=numpy.ones(shape=(len(numpy.load("carb_RNG200my_t.npy")),3))
carb_filt[:,0]=numpy.load("carb_RNG200my_t.npy") #time values for carbonates
carb_filt[:,1]=numpy.load("carb_RNG200my_v.npy") #isotope values for carbonates
carb_filt[:,2]=numpy.load("carb_RNG200my_var.npy") #isotope values for carbonates

obs_carb=carb_filt[:,1]
time_carb=carb_filt[:,0]
er_carb=carb_filt[:,2]

## Load binned organic carbon isotope data
org_filt=numpy.ones(shape=(len(numpy.load("org_200my_t.npy")),3))
org_filt[:,0]=numpy.load("org_200my_t.npy") #time values for organics
org_filt[:,1]=numpy.load("org_200my_v.npy") #isotope values for organics
org_filt[:,2]=numpy.load("org_200my_var.npy") #isotope values for organics

org_filt=numpy.ones(shape=(len(numpy.load("org_RNG200my_t.npy")),3))
org_filt[:,0]=numpy.load("org_RNG200my_t.npy") #time values for organics
org_filt[:,1]=numpy.load("org_RNG200my_v.npy") #isotope values for organics
org_filt[:,2]=numpy.load("org_RNG200my_var.npy") #isotope values for organics

obs_org=org_filt[:,1]
time_org=org_filt[:,0]
er_org=org_filt[:,2]
#########################################


global Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,cons_dict

# Define constants
n = 1.7 #exponent for carbonate precipitation
F_out_mantle_mod=3e12 #Modern mantle outgassing of carbon, mol/yr
new_add_Ca=0.0    #Redundant
F_meta_mod=3e12 #Modern metamorphic outgassing of carbon, mol/yr
del_carb0=-5.5  #Initial carbonate isotopic ratio (4.1 Ga), per mil
del_org0=-5.5   #Initial organic isotopic ratio (4.1 Ga), per mil
del_mantle0=-5.5 #Initial mantle isotopic ratio (4.1 Ga), per mil
F_sub_mod=6e12 #Modern subduction flux of carbon, mol/yr
e2=0.5 #Organic subduction efficiency
Bio_f = 0.99999 #Redundant
coef_for_diss = 0.25 #Seafloor weathering dissolution coefficient

## Define likelihood function for inverse analysis
def LnLike(x):
   
    ## Define prior ranges for unknown parameters:
    j1=x[0] # Archean organic burial fraction
    if (x[0] > 0.5 ) or (x[0]<0.01): 
        print ("oops, didn't work 0")
        return -numpy.inf    
    
    j2=x[1] # Multiplicative change organic burial Archean -> Proterozoic
    if (x[1] > 5.0 ) or (x[1]<0.5): 
        print ("oops, didn't work 1")
        return -numpy.inf  
        
    j3=x[2] # Multiplicative change organic burial Proterozoic -> Phanerozoic
    if (x[2] > 5.0 ) or (x[2]<0.5): 
        print ("oops, didn't work 2")
        return -numpy.inf

    if j1*j2*j3 >1.0: # Check that modern forg does not exceed 1:
        print (">1 forg")
        return -numpy.inf
    
    org_weath_mod0=x[3] #Modern organic weathering flux, mol C/yr
    #if (x[3] > 8e12 ) or (x[3]< 2e12): # Sensitivity test
    if (x[3] > 5e12 ) or (x[3]< 2e12): # Nominal range
        print ("oops, didn't work 3")
        return -numpy.inf 
    
    thermo_weath=x[4] #Modern thermogenic flux, mol C/yr
    if (x[4] > 4e12 ) or (x[4]< 1e12): # Nominal range
        print ("oops, didn't work 4")
        return -numpy.inf            
    
    protero_O2=10**x[5] #Proterozoic oxygen, PAL
    if (x[5] > -1 ) or (x[5]< -3 ): 
        print ("oops, didn't work 5")
        return -numpy.inf
    
    ppCO2_00 = 10**x[6] #Initial (late Hadean) CO2, bar
    if (x[6] > 1.5) or (x[6]<-2.5):
        print ("oops didin't work 6")
        return -numpy.inf

    climp = x[7] # CO2 silicate weathering feedback exponent
    if (x[7] > 0.5) or (x[7]<0.1):
        print ("oops didin't work 7")
        return -numpy.inf

    tdep_weath = x[8] # T-dependence of silicate weathering
    if (x[8] > 40.0) or (x[8]< 10.0):
        print ("oops didin't work 8")
        return -numpy.inf

    beta = x[9] #Exponent controlling dependence of outgassing on heatflow
    if (x[9] > 2.0) or (x[9]< 0.0):
        print ("oops didin't work 9")
        return -numpy.inf

    mm = x[10] #Exponent controlling dependence of spreading rate on heatflow
    if (x[10] > 2.0) or (x[10] < 1.0e-3):
        print ("oops didin't work 10")
        return -numpy.inf

    n_out = x[11] #Exponent controlling heatflow evolution
    if (x[11] > 2.0) or (x[11] < 0.0):
        print ("oops didin't work 11")
        return -numpy.inf

    Ebas = x[12] #Effective activation energy seafloor weathering
    if (x[12] > 100000) or (x[12] < 60000.0):
        print ("oops didin't work 12")
        return -numpy.inf

    multiplicative = x[13] #Change in reducing power volcanic gases
    #if (x[13] > 4.0) or (x[13] < 0.0): #Reduced mantle sensitivity test
    if (x[13] > 0.000004) or (x[13] < 0.0): ## Nominal model
        print ("oops didin't work 13")
        return -numpy.inf

    F_carbw = x[14] #Modern carbonate weathering flux, mol/yr
    #if (x[14] > 15e12) or (x[14] < 7e12): #Sensitivity test
    if (x[14] > 25e12) or (x[14] < 7e12): #Nominal model
        print ("oops didin't work 14")
        return -numpy.inf
   
    e1 = x[15] # Modern carbonate subduction efficiency
    if (x[15] > 0.8) or (x[15] < 0.2):
        print ("oops didin't work 15")
        return -numpy.inf

    lfrac = x[16] # Land fraction parameter
    #if (x[16] >  2.0) or (x[16] < 1.1111): #Sensitivity test
    #if (x[16] >  0.8) or (x[16] < 0.4): #Sensitivity test
    if (x[16] >  2.0) or (x[16] < 0.4): #Nominal 
        print ("oops didin't work 16")
        return -numpy.inf
    
    growth_timing = x[17] # Timing of growth of continents, Ga
    if (x[17] > 3.0) or (x[17] < 2.0):
        print ("oops didin't work 17")
        return -numpy.inf
    
    pH_00 = x[18] # Initial (late Hadean) ocean pH
    if (x[18] > 8.5) or (x[18] < 5.0):
        print ("oops didin't work 18")
        return -numpy.inf

    R_carb_00 = (10e21)*(10**x[19]) # Initial (late Hadean) crustal carbonate reservoir, mol C
    if (x[19] > 0.0) or (x[19] < -4.0):
        print ("oops didin't work 19")
        return -numpy.inf

    R_mantle_00 = x[20] # Initial (late Hadean) mantle carbon reservoir, mol C
    if (x[20] > 4e22) or (x[20] < 0.1e22):
        print ("oops didin't work 20")
        return -numpy.inf

    Mantle_temp_mod = x[21] # Initial mantle temperature (K)
    if (x[21] > 1800) or (x[21] < 1300):
        print ("oops didin't work 21")
        return -numpy.inf
    Ccrit0 = (3.125e-3*80.0+835.5 - 285.0)/(Mantle_temp_mod-285.0)
    carb_eff_modifier0 = (Ccrit0 - 0.3)/(0.6-0.3) #carbonate subduction efficiency 
    if carb_eff_modifier0 > 1-e1: ## Check to ensure carbonate subduction more efficient than organic subduction
        print ("carbonates better subducted than organics!")
        return -numpy.inf        

    ## Attempt to run forward model
    try:    
        [[out0,out1,out2,out3,time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array,volc_array,F_meta_ar,F_sub_ar,F_out_mantle_ar,out4,out5,out6,F_meta_carb_ar,F_sub_carb_ar,F_meta_org_ar,F_sub_org_ar,org_weath_ar,F_orgB_ar,out7,out8,out9,out10,del_crust_ar,del_out_ar,del_inputs_ar,carbw_ar,del_org_buried_ar,del_carb_precip_ar,pO2_ar,reduced_early_ar,F_sub_arc_ar,del_arc_volc_ar,del_meta_volc_ar,del_sub_ar,one_minus_carb_eff_modifier_ar,Kox_array],imbalance]=Precambrian(Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,cons_dict,multiplicative,Mantle_temp_mod)																										
                                                                                                                                                                                 
        print("success")
    except:
        print ("oops, didn't work")
        return -numpy.inf

    if numpy.min(Fs_array)< 0:
        print ("Negative silicate weathering")
        return -numpy.inf # if negative silicate weathering discard

    # Discard non-physical, or nan solutions:
    if (numpy.max(out0)>1e10)or(numpy.min(out0)<=0):
        print ("Kill because crazy high 0")
        return -numpy.inf #added
    if (numpy.max(out1)>1e10)or(numpy.min(out1)<=0):
        print ("Kill because crazy high 1")
        return -numpy.inf     
    if (numpy.max(out2)>1e10)or(numpy.min(out2)<=0):
        print ("Kill because crazy high 2")
        return -numpy.inf #added
    if (numpy.max(out3)>1e10)or(numpy.min(out3)<=0):
        print ("Kill because crazy high 3")
        return -numpy.inf      
    if numpy.any(numpy.isnan([out0,out1,out2,out3,time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array,volc_array,F_meta_ar,F_sub_ar,F_out_mantle_ar,out4,out5,out6,F_meta_carb_ar,F_sub_carb_ar,F_meta_org_ar,F_sub_org_ar,org_weath_ar,F_orgB_ar,out7,out8,out9,out10,del_crust_ar,del_out_ar,del_inputs_ar,carbw_ar,del_org_buried_ar,del_carb_precip_ar,pO2_ar])):
        print ("Kill because nan")
        return -numpy.inf

    mass_conserve = out4+out6+out5+out0*1.35e21
    print (numpy.max(mass_conserve)-numpy.min(mass_conserve),"mass conserve")
    balance_masses=numpy.max(mass_conserve)-numpy.min(mass_conserve)
    if (balance_masses>1e12): #Discard if didn't conserve mass
        print ("Kill because no mass conserve")
        return -numpy.inf


    ## Interpolated carbonate and organic isotope values for comparison to data
    carb_interp = interp1d(time, del_carb_precip_ar) 
    org_interp = interp1d(time, del_org_buried_ar) 
    model_carb=carb_interp(-time_carb*1e6)
    model_org=org_interp(-time_org*1e6)

    ## Get observed values for model parameters to compute likelihood
    modern_CO2 = numpy.log10(CO2_array_o[len(CO2_array_o)-1]) #model modern pCO2
    true_mod_CO2 = cons_dict['true_mod_CO2'] #observed modern pCO2
    er_CO2 = cons_dict['er_CO2']  # error in pCO2
    modern_Temp = Tsurf_array[len(CO2_array_o)-1] #model modern Temp
    true_mod_Temp = cons_dict['true_mod_Temp'] #observed modern Temp
    er_modeTemp = cons_dict['er_modeTemp'] # error in Temp
    model_modern_forg = j1*j2*j3 #model modern forg
    true_mod_forg = cons_dict['true_mod_forg'] #observed modern forg
    er_true_mod_forg = cons_dict['er_true_mod_forg'] #error forg
    true_mod_pH = cons_dict['true_mod_pH'] #observed modern pH
    er_true_mod_pH = cons_dict['er_true_mod_pH'] #error modern pH
    model_modern_pH = pH_array_o[999] #model modern pH

    # Avice et al. error in mantle outgassing at 3.3 Ga
    Avice_outgassing = cons_dict['Avice_outgassing'] # Observed outgassing
    error_Avice_outgassing = cons_dict['error_Avice_outgassing'] # error in observed outgassing
    actual_3_3 = F_out_mantle_ar[194]/F_out_mantle_ar[999] # model outgassing
    
    Archean_mantle_del13C = cons_dict['Archean_mantle_d13C'] #Observed Archean mantle d13C
    error_Archean_mantle_del13C = cons_dict['er_Archean_mantle_d13C'] #Error Archean mantle d13C
    model_Archean_d13Cmantle = numpy.median(out9[49:365]) #Model Archean mantle d13C
    
    Modern_mantle_del13C = cons_dict['Archean_mantle_d13C'] #Observed modern mantle d13C
    error_Modern_mantle_del13C = cons_dict['er_Archean_mantle_d13C']/9.0 #Error modern mantle d13C
    model_Modern_d13Cmantle = numpy.median(out9[999]) #Model modern mantle d13C

    weath_burial_const = 0
    counter_ratio = 0
    for ii in range(0,len(org_weath_ar)): ## Add constraint to ensure organic burial exceeds organic weathering
        if (org_weath_ar[ii] > F_orgB_ar[ii])and(F_orgB_ar[ii]>0):
            counter_ratio = counter_ratio + 1
            weath_burial_const = weath_burial_const + ((org_weath_ar[ii] - F_orgB_ar[ii])/F_orgB_ar[ii])**2
    weath_burial_const = 100 * weath_burial_const

    Kox = numpy.median(Kox_array[49:365]) #Model Kox
    er_Kox = 0.20 #Error in Kox
    if Kox < 1: # Constrain Kox to be less than 1, within error
        Kox_term = 0.0
    else:
        Kox_term = (Kox - 1)**2/er_Kox**2
 
    ## Define observed modern reservoirs and their uncertainties 
    mantle_mod_true= cons_dict['mantle_mod_true'] # mantle C
    er_mantle_true = cons_dict['er_mantle_true'] 
    Rcrust_mod_true= cons_dict['Rcrust_mod_true'] # crustal carbonates
    er_Rcrust_true = cons_dict['er_Rcrust_true']
    Rorg_mod_true= cons_dict['Rorg_mod_true'] # crustal organics
    er_Rorg_true = cons_dict['er_Rorg_true']

    Res_model_mantle = out5[999] #Model modern mantle reservoir
    Res_model_carb = out4[999] #Model modern carbonate reservoir
    Res_model_org = out6[999] #Model modern organic reservoir
    weight_c = 1.0 #Weighting (redundant)
    weight_o = 1.0 #Weighting (redundant)

    ## Now calculate log-likelihood 
    # Constant log terms:   
    log_terms= weight_o*numpy.sum(numpy.log(2*numpy.pi*er_org**2)) + weight_c*numpy.sum(numpy.log(2*numpy.pi*er_carb**2))+ numpy.log(2*numpy.pi*er_CO2**2) + numpy.log(2*numpy.pi*er_modeTemp**2)+ 0*numpy.log(2*numpy.pi*er_true_mod_forg **2) + numpy.log(2*numpy.pi*error_Avice_outgassing **2) + numpy.log(2*numpy.pi*er_mantle_true**2) + numpy.log(2*numpy.pi*er_Rcrust_true**2) + numpy.log(2*numpy.pi*er_Rorg_true**2) + numpy.log(2*numpy.pi*er_Kox**2) + numpy.log(2*numpy.pi*er_true_mod_pH**2) + numpy.log(2*numpy.pi*error_Archean_mantle_del13C**2) + numpy.log(2*numpy.pi*error_Modern_mantle_del13C**2)

    # Error terms:
    log_like =  -0.5 * (  weight_c*numpy.sum((obs_carb-model_carb)**2/er_carb**2) +  weight_o* numpy.sum((obs_org-model_org)**2/er_org**2) +  numpy.sum((numpy.log10(true_mod_CO2)-modern_CO2)**2/er_CO2**2) +  numpy.sum((true_mod_Temp- modern_Temp)**2/er_modeTemp**2)  +  0*numpy.sum((model_modern_forg - true_mod_forg)**2/er_true_mod_forg **2) +numpy.sum((Avice_outgassing - actual_3_3)**2/error_Avice_outgassing**2) +(Res_model_mantle - mantle_mod_true)**2/er_mantle_true**2 + (Res_model_carb - Rcrust_mod_true)**2/er_Rcrust_true**2  + (Res_model_org - Rorg_mod_true)**2/er_Rorg_true**2 + Kox_term + (true_mod_pH - model_modern_pH)**2/er_true_mod_pH**2 + (Archean_mantle_del13C - model_Archean_d13Cmantle)**2/error_Archean_mantle_del13C**2 + (Modern_mantle_del13C - model_Modern_d13Cmantle)**2/error_Modern_mantle_del13C**2 + weath_burial_const + log_terms)

    if numpy.isnan(log_like):
        print("log like is NAN")
        return -numpy.inf
    
    return log_like


f = open("../constraint_dict_revised.txt",mode='r',encoding = 'utf-8-sig') #Contains observed/true values
cons_dict = eval(f.read())
nwalk = 500 #Number of walkers in MCMC calculation
ndim = 22 #Number of unknown paramteres
nsteps = 10000 #Number of steps per walker, total forward model calls = nwalk*nsteps
number_of_cores = 50 #Threads used to parallelize calculation

## Initialize walkers:
p0 = numpy.vstack([[0.01+0.49*numpy.random.random() for i in range(nwalk)],
                [0.5+4.5*numpy.random.random() for i in range(nwalk)],
                [0.5+4.5*numpy.random.random() for i in range(nwalk)],
                [2e12+3e12*numpy.random.random() for i in range(nwalk)],
                [1e12+3e12*numpy.random.random() for i in range(nwalk)],
                [-3.0+2.0*numpy.random.random() for i in range(nwalk)],
                [-2.5+4.0*numpy.random.random() for i in range(nwalk)],
                [0.1+0.4*numpy.random.random() for i in range(nwalk)],
                [10.0+30.0*numpy.random.random() for i in range(nwalk)],
                [0.0+2.0*numpy.random.random() for i in range(nwalk)],
                [0.0+2.0*numpy.random.random() for i in range(nwalk)],
                [0.0+2.0*numpy.random.random() for i in range(nwalk)],
                [60000+40000*numpy.random.random() for i in range(nwalk)],
                [0.000004*numpy.random.random() for i in range(nwalk)],
                #[7e12+8e12*numpy.random.random() for i in range(nwalk)],
                #[10e12+15e12*numpy.random.random() for i in range(nwalk)],
                [7e12+18e12*numpy.random.random() for i in range(nwalk)],
                [0.2+0.6*numpy.random.random() for i in range(nwalk)],
                #[0.4+0.4*numpy.random.random() for i in range(nwalk)],
                #[1.11111+0.888888*numpy.random.random() for i in range(nwalk)],
                [0.4+1.6*numpy.random.random() for i in range(nwalk)],
                [2.0+1.0*numpy.random.random() for i in range(nwalk)],
                [5.0+3.5*numpy.random.random() for i in range(nwalk)],
                [-4.0 + 4.0*numpy.random.random() for i in range(nwalk)],
                [0.1e22 + 3.9e22*numpy.random.random() for i in range(nwalk)],
                [1300 + 500*numpy.random.random() for i in range(nwalk)]]).T


## Run emcee sampler:
sampler = emcee.EnsembleSampler(nwalk, ndim, LnLike,threads=number_of_cores) #was before 14, 8 safer maybe
                                #args = (t, y, error),
                                #pool = emcee.interruptible_pool.InterruptiblePool()) #for parallel
pos, lnprob, rstate=sampler.run_mcmc(p0, nsteps)


## Save chains
numpy.save('chain_output',sampler.chain[:,:,:]) 
numpy.save('lnprob_output',sampler.lnprobability)

# Load chains (if not re-running entire code, start here)
chain=numpy.load('chain_output.npy') 
lnprob=numpy.load('lnprob_output.npy')  

# Plot the chains
fig, ax = pl.subplots(15) 
for nnn in range(nwalk):
  ax[0].plot(chain[nnn,:,0])
  ax[1].plot(chain[nnn,:,1])
  ax[2].plot(chain[nnn,:,2]) 
  ax[3].plot(chain[nnn,:,3]) 
  ax[4].plot(chain[nnn,:,4]) 
  ax[5].plot(chain[nnn,:,5]) 
  ax[6].plot(chain[nnn,:,6]) 
  ax[7].plot(chain[nnn,:,7]) 
  ax[8].plot(chain[nnn,:,8]) 
  ax[9].plot(chain[nnn,:,9]) 
  ax[10].plot(chain[nnn,:,10]) 
  ax[11].plot(chain[nnn,:,11]) 
  ax[12].plot(chain[nnn,:,12]) 
  ax[13].plot(chain[nnn,:,13])


fig, ax = pl.subplots(12) 
for nnn in range(nwalk):
    if (numpy.isinf(lnprob[nnn,nsteps-1])):
        dummy = 45
    else:
        ax[0].plot(chain[nnn,:,0])
        ax[1].plot(chain[nnn,:,1])
        ax[2].plot(chain[nnn,:,2])
        ax[3].plot(chain[nnn,:,3])
        ax[4].plot(chain[nnn,:,4]) 
        ax[5].plot(chain[nnn,:,5]) 
        ax[6].plot(chain[nnn,:,6]) 
        ax[7].plot(chain[nnn,:,7]) 
        ax[8].plot(chain[nnn,:,8]) 
        ax[9].plot(chain[nnn,:,9]) 
        ax[10].plot(chain[nnn,:,10]) 
        ax[11].plot(chain[nnn,:,11]) 

#find highest likelihood run
logprob=numpy.array(lnprob)
values=chain
ii,jj = numpy.unravel_index(logprob.argmax(), logprob.shape)
print ("indeces for best",ii,jj)
print ("loglikelihood and values",logprob[ii,jj],values[ii,jj,:])
print ("Check",LnLike(values[ii,jj,:]))

# Create corner plots
production = chain[:,:,:]
#production = chain[:,2000:,:] #discarding the first 2000 steps as burn-in
#timmed_lnprob = lnprob[:,2000:] #discarding the first 2000 steps as burn-in

production = chain[:,:,:] #no discarding burn-in
timmed_lnprob = lnprob[:,:] #no discarding burn-in

s = production.shape
ss= timmed_lnprob.shape
flatchain = production.reshape(s[0] * s[1], s[2])
flatlnprob= timmed_lnprob.reshape(s[0] * s[1])

flatchain = flatchain[numpy.logical_not(numpy.isinf(flatlnprob))]
flatlnprob = flatlnprob[numpy.logical_not(numpy.isinf(flatlnprob))]

# Corner plots
copy_flatchain = numpy.copy(flatchain)
copy_flatchain[:,3] = copy_flatchain[:,3]/1e12 #Convert Tmol/yr
copy_flatchain[:,4] = copy_flatchain[:,4]/1e12 #Convert Tmol/yr
copy_flatchain[:,20] = copy_flatchain[:,20]/1e21 #Convert Tmol/yr

titles=numpy.array(['$f_{org}^{Archean}$ = $j_1$', '$f_{org}^{ Protero}$/$f_{org}^{ Archean}$ = $j_2$', '$f_{org}^{ Phanero}$/$f_{org}^{ Protero}$ = $j_3$', '$F_{oxid}^{mod}$', '$F_{thermo}^{mod}$','$log(pO_2^{Protero}/pO_2^{mod})$', '$log(pCO_2^{init})$', '$n_{out}$', '$R_{mantle}^{init}$'])
corner.corner(copy_flatchain[:,[0,1,2,3,4,5,6,11,20]], quantiles=[0.16,0.5,0.84], labels=titles, show_titles='True', plot_contours='True')

titles=numpy.array(['$f_{org}^{Archean}$ = $j_1$', '$f_{org}^{ Protero}$/$f_{org}^{ Archean}$ = $j_2$', '$f_{org}^{ Phanero}$/$f_{org}^{ Protero}$ = $j_3$', 'MantleT', 'e1','Multiplicative', 'Fcarbw', 'OrgWeath', 'ThermoWeath'])
corner.corner(copy_flatchain[:,[0,1,2,21,15,13,14,3,4]], quantiles=[0.16,0.5,0.84], labels=titles, show_titles='True', plot_contours='True')

pl.figure()
pl.hist(flatchain[:,1]*flatchain[:,2],bins=30,normed=True)
pl.xlabel('Relative change in forg over Earth history')
pl.ylabel('Probability density')

ab, bc, cd,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,rs,st,tu,uv,vw = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*numpy.percentile(flatchain, [16, 50, 84],axis=0)))
print ("median values with errors", numpy.array([ab, bc, cd,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,rs,st,tu,uv,vw]))
print ("confidence intervals")
map(lambda v: (v[0], v[1], v[2]),zip(*numpy.percentile(flatchain, [5, 50, 95],axis=0)))

## Create time-evolution plots (will need to re-run forward model on represenative sample)
from plotting_everything import dist_plotter #load plotting script

mega_output=[]
spread_output=[]
carbw_factor=[]
sil_change=[]
seafloor_change=[]

inf_counter=0
weird_counter=0
total_sampled=1000 # Run forward model for 1000 parameter values
for i in range(0,total_sampled):
    j=numpy.random.randint(len(flatchain))
    x_ex=flatchain[j]
    ln_like = flatlnprob[j]
    print (ln_like,j)
    if numpy.isinf(ln_like):
        inf_counter=inf_counter+1
        print ("Whoops, -inf likelihood")
    else: 
        print ("seems to be OK")
        try:
            [outputs,imblanace]=Precambrian(Bio_f,10**x_ex[6],n,x_ex[7],x_ex[8],F_meta_mod,F_out_mantle_mod,x_ex[18],x_ex[16],(10e21)*(10**x_ex[19]) , del_carb0,x_ex[14], del_mantle0, x_ex[20], F_sub_mod, coef_for_diss,x_ex[9], x_ex[11],x_ex[10],x_ex[17],new_add_Ca,x_ex[12], x_ex[15] ,10**x_ex[5], del_org0, e2, x_ex[0], x_ex[1], x_ex[2], x_ex[3], x_ex[4], cons_dict,x_ex[13],x_ex[21])
            if numpy.any(numpy.isnan(outputs)):
                print ("Kill because nan")
                weird_counter=weird_counter+1
            else:
                spread_output.append( 2)
                mega_output.append(outputs)
        except:
            print ("something weird happening")
            weird_counter=weird_counter+1

print(numpy.shape(chain))
print(numpy.shape(flatchain))
print( "Acceptance fraction",sampler.acceptance_fraction)
print (inf_counter," of ",total_sampled)
print ("weird counter",weird_counter)
    
## Plot these 1000 outputs:
mega_output=numpy.array(mega_output)
spread_output=numpy.array(spread_output)
dist_plotter(mega_output,spread_output,"y")





