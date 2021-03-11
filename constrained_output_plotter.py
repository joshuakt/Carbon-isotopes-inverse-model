import emcee
import numpy
import matplotlib.pyplot as pl
from corner import corner
import corner


#from reverse_functions_no_pore_modified_new import Precambrian
from reverse_functions_no_pore_modified_new_subduct import Precambrian
from scipy.interpolate import interp1d

###########################

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
#er_carb[-1] = er_carb[-1] / 50
#er_carb[-2] = er_carb[-2] / 50
#######################

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

##############################

#preinudsmod=280.0
#time_CO2=numpy.array([65,75,85,95])
#obs_CO2=numpy.array([441,1034,1220,1528])/preinudsmod
#er_CO2=numpy.array([184,535,428,95])/preinudsmod #these currently running

########################

#obs_prec_zero=2.35e12 #need spread modifier later
#time_prec=99
#er_prec_zero=0.75e12 # need spread modifier late

########################

#obs_omega_o=2.48#4.1 #from figure
#time_omega=99
#er_omega_o=0.3#1.3 #from nowhere in particular

#########################

#obs_pH=7.72 #average of middle values
#time_pH=50
#er_pH=0.1 # ~std err

##########################
global Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,cons_dict

#ppCO2_00=1.2 # bar
#n=1.7 #numpy.random.uniform(1.0,2.5)
n = 1.7
F_out_mantle_mod=3e12 #numpy.random.uniform(0.5,1.5) #0.7 #IRRELEVANT CAN REPLACE
#tdep_weath=20 #numpy.random.uniform(10.,40.)
#climp=0.2# numpy.random.uniform(0.2,0.5) #0.2-0.5
#pH_00=6.2
#lfrac=1.5 #numpy.random.uniform(1.111,2)#from 1.111-10% to 2-50% 4-75% Archean 1.1111 2.0
#growth_timing=2.5 #numpy.random.uniform(2.0,3.0)
new_add_Ca=0.0    #numpy.random.uniform(0,160.0) #IRRELEVANT CAN REPLACE
F_meta_mod=3e12 #numpy.random.uniform(0.25e12,0.65e12) ## CURRENTLY INTIAL NOT MODERN BECAUSE DIFFICULT WITH pH_p
del_carb0=-5.5  #2.3 numpy.random.uniform(1.0,2.0) 
del_org0=-5.5   #30 #numpy.random.uniform(-24,-30)
del_mantle0=-5.5 #numpy.random.uniform(-5,-6)
#F_carbw=10e12 #numpy.random.uniform(10e12,14e12) #see rabbit hole email CAN BE EITHER MODERN OR INITIAL COND.
#R_mantle_00=100000e18 #numpy.random.uniform(50000e18 ,150000e18 )
#R_carb_00=0.008*7000e18# numpy.random.uniform(1000e18,10000e18)#
F_sub_mod=6e12 #numpy.random.uniform(0.8,1.4)
#coef_for_diss=0.25 #numpy.random.uniform(0.0,0.5)
#beta=1.5 #numpy.random.uniform(0.0,2.0)
#mm=1.0#2.0 #numpy.random.uniform(1.0,2.0)
#n_out=0.66#.73#.5# numpy.random.uniform(0.0,0.73)#0.73#.4 #numpy.random.uniform(0.0,0.73) #for heatflow!
#Ebas=90000 #numpy.random.uniform(60000.,100000.)
#e1=0.3#0.5#0.3 #organic recycle to arc
e2=0.5#0.5#0.7 # carbonate recycle to arc
#protero_O2=0.01 #10**(numpy.random.uniform(-4,-1)) #PAL
Bio_f = 0.99999#0.5
coef_for_diss = 0.25


# Randomize some data
#t = np.linspace(0, 10, 100)
#m = 5.
#b = 3.
#y_true = m * t + b
#y = y_true + 2. * np.random.randn(100)
#error = np.ones(100) * 2.

f = open("../constraint_dict_revised.txt",mode='r',encoding = 'utf-8-sig')
cons_dict = eval(f.read())
# Set up the sampler
nwalk = 500 #0 # #100 x 1000 for <2 hours laptop 600x1000 drafts, no burnin, usually 1000
ndim = 22 #changed from 16 befor PG nonsesne
nsteps = 20000 #00 #10000

#nwalk = 300 #0 # #100 x 1000 for <2 hours laptop 600x1000 drafts, no burnin, usually 1000
#ndim = 22 #changed from 16 befor PG nonsesne
#nsteps = 6000 #00 #10000



chain=numpy.load('chain_TEST22b.npy') #TEST is for 100 timesteps 100x1000, 5 var
lnprob=numpy.load('lnprob_TEST22b.npy')  #TEST3 for 6 var, 1000x1000




#find highest likelihood run
logprob=numpy.array(lnprob)
values=chain

#sm=numpy.load('sampler.npy')

# Plot the corner plot, discarding the first 100 steps as burn-in
production = chain[:,:,:]
production = chain[:,2000:,:]
timmed_lnprob = lnprob[:,2000:]

#production = chain[:,:,:]
#timmed_lnprob = lnprob[:,:]

s = production.shape
ss= timmed_lnprob.shape
flatchain = production.reshape(s[0] * s[1], s[2])
flatlnprob= timmed_lnprob.reshape(s[0] * s[1])
#flatchain=sampler.flatchain
#corner.corner(flatchain, quantiles=[0.16, 0.5, 0.84],labels=["CO2-dep", "T-dep", "Weatherability","climate","outgas_change","CWF","initial_out","initial_cweath","circulation_time","n_precip","alt_frac","deep_grad","pHdep_sea","Eact","fpel","beta"],truths=values[ii,jj,:])


## Throw away large modern organic weathering fluxes
#for k in range(0,len(flatlnprob)):
#    if flatchain[k,3] > 5e12:
#        flatlnprob[k] = numpy.inf


flatchain = flatchain[numpy.logical_not(numpy.isinf(flatlnprob))]
flatlnprob = flatlnprob[numpy.logical_not(numpy.isinf(flatlnprob))]

# more limited corner plots
#corner.corner(flatchain[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]], quantiles=[0.16, 0.5, 0.84],labels=["forg init.", "GOE increase", "NOE increase","O2 org. weath.","Thermo. org. weath.","Protero. O2","pCO2 init","climp","Te_weath","beta","mm","n_out","Ebas","coef sea","Fcarbw"])#,truths=values[ii,jj,:])
#pl.figure(figsize=(40,30))

#corner.corner(flatchain[:,[11,12,13,14,15,16,17,18,19,20]], quantiles=[0.16, 0.5, 0.84],labels=["n_out","Ebas","coef sea","Fcarbw","n","lfrac","growth time.","pH_00","Rcarb_00","Rmantly_00"])#,truths=values[ii,jj,:])
#pl.figure(figsize=(40,20))

copy_flatchain = numpy.copy(flatchain)
copy_flatchain[:,3] = copy_flatchain[:,3]/1e12
copy_flatchain[:,4] = copy_flatchain[:,4]/1e12
copy_flatchain[:,20] = copy_flatchain[:,20]/1e21
titles=numpy.array(['$f_{org}^{Archean}$ = $j_1$', '$f_{org}^{ Protero}$/$f_{org}^{ Archean}$ = $j_2$', '$f_{org}^{ Phanero}$/$f_{org}^{ Protero}$ = $j_3$', '$F_{oxid}^{mod}$', '$F_{thermo}^{mod}$','$log(pO_2^{Protero}/pO_2^{mod})$', '$log(pCO_2^{init})$', '$n_{out}$', '$R_{mantle}^{init}$'])
corner.corner(copy_flatchain[:,[0,1,2,3,4,5,6,11,20]], quantiles=[0.16,0.5,0.84], labels=titles, show_titles='True', plot_contours='True')#, truths=trutharr)
#corner.corner(flatchain[:,[0,1,2,3,4,5,6,11,20]], quantiles=[0.16, 0.5, 0.84],labels=["forg init.", "GOE increase", "NOE increase","O2 org. weath.","Thermo. org. weath.","Protero. O2","pCO2 init","n_out","Rmantly_00"])#,truths=values[ii,jj,:])

titles=numpy.array(['$f_{org}^{Archean}$ = $j_1$', '$f_{org}^{ Protero}$/$f_{org}^{ Archean}$ = $j_2$', '$f_{org}^{ Phanero}$/$f_{org}^{ Protero}$ = $j_3$', 'MantleT', 'e1','Multiplicative', 'Fcarbw', 'OrgWeath', 'ThermoWeath'])
corner.corner(copy_flatchain[:,[0,1,2,21,15,13,14,3,4]], quantiles=[0.16,0.5,0.84], labels=titles, show_titles='True', plot_contours='True')

import pylab
pylab.show()
import pdb
pdb.set_trace()


## Show the plots
#from plotting_everything_no_pore_new import mc_plotter_spread,dist_plotter
from plotting_everything_no_pore_final import mc_plotter_spread,dist_plotter
#[array_outputs,imbalance]=Precambrian(W,10**values[ii,jj,6],n,values[ii,jj,7],values[ii,jj,8],mod_sea,alt_frac,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,deep_grad,coef_for_diss,values[ii,jj,9],values[ii,jj,11],values[ii,jj,10],growth_timing,new_add_Ca,Ebas,e1,10**values[ii,jj,5],del_org0,e2,values[ii,jj,0],values[ii,jj,1],values[ii,jj,2],values[ii,jj,3],values[ii,jj,4])

## commented out
##spread_best= ( (values[ii,jj,6]+values[ii,jj,4]*values[ii,jj,6])/values[ii,jj,6] )**values[ii,jj,15]
#spread_best= 2
#mc_plotter_spread(array_outputs,"n",0,spread_best)
#
#import pylab
#pylab.figure(figsize=(30,15))
#legend_counter=0
#for x_ex in flatchain[numpy.random.randint(len(flatchain), size=100)]:
#    print x_ex
#    [outputs,imblanace]=Precambrian(W,ppCO2_00,n,climp,tdep_weath,mod_sea,alt_frac,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,deep_grad,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,10**x_ex[5],del_org0,e2,x_ex[0],x_ex[1],x_ex[2],x_ex[3],x_ex[4])
#    #mega_output.append(outputs)
#    sp=4
#    mc_plotter_spread(outputs,"y",legend_counter,sp)
#    legend_counter=legend_counter+1
## end commented out

mega_output=[]
spread_output=[]
carbw_factor=[]
sil_change=[]
seafloor_change=[]

import pdb
inf_counter=0
weird_counter=0
total_sampled=1000
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
            elif outputs[38][999] < -7.0:
                weird_counter=weird_counter+1
            else:
                spread_output.append( 2)
                mega_output.append(outputs)
        except:
            print ("something weird happening")
            weird_counter=weird_counter+1

print(numpy.shape(chain))
print(numpy.shape(flatchain))
print (inf_counter," of ",total_sampled)
print ("weird counter",weird_counter)
    
mega_output=numpy.array(mega_output)
spread_output=numpy.array(spread_output)
dist_plotter(mega_output,spread_output,"y")





