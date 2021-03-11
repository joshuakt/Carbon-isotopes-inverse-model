#! /usr/local/bin/python
import numpy
import pylab

import scipy.stats
from scipy.interpolate import interp1d

def dist_plotter(all_output,spread_output,sd):   

    f = open("../constraint_dict_revised.txt",mode='r',encoding = 'utf-8-sig')
    cons_dict = eval(f.read())

    confidence_pH_o=scipy.stats.scoreatpercentile(all_output[:,5,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_pH_p=scipy.stats.scoreatpercentile(all_output[:,7,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)

    pylab.figure(figsize=(30,15))
    pylab.subplot(3, 4, 1)
    true_mod_pH = cons_dict['true_mod_pH']
    er_true_mod_pH = cons_dict['er_true_mod_pH']
    pylab.plot(all_output[0,4,:],confidence_pH_o[1],'r',label='ocean')
    pylab.plot(all_output[0,4,:],confidence_pH_p[1],'b',label='pore space')
    pylab.fill_between(all_output[0,4,:], confidence_pH_o[0], confidence_pH_o[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_pH_p[0], confidence_pH_p[2], color='blue', alpha='0.4')
    pylab.xlabel('Time (yr)')
    pylab.ylabel('pH')

    ############################
    import matplotlib
    ppCO2=10**-6
    confidence_CO2o=scipy.stats.scoreatpercentile(all_output[:,6,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.subplot(3, 4, 2)
    pylab.plot(all_output[0,4,:],confidence_CO2o[1]/ppCO2,'r',label='RCO2')
    pylab.fill_between(all_output[0,4,:], confidence_CO2o[0]/ppCO2, confidence_CO2o[2]/ppCO2, color='red', alpha='0.4')
    true_mod_CO2 = cons_dict['true_mod_CO2']/ppCO2
    er_CO2 = cons_dict['er_CO2']
    positive_er = ( -10**( numpy.log10(cons_dict['true_mod_CO2']) - er_CO2 ) + cons_dict['true_mod_CO2'] ) /ppCO2
    negative_er = ( 10**( numpy.log10(cons_dict['true_mod_CO2']) + er_CO2 ) - cons_dict['true_mod_CO2'] ) /ppCO2
    matplotlib.pyplot.errorbar(all_output[0,4,999],true_mod_CO2,yerr=numpy.transpose([[negative_er,positive_er]]),color='k',marker='o')
    pylab.xlabel('Time (yr)')
    pylab.ylabel('CO2 relative to modern')
    pylab.legend(loc=2)
    
    #########################
    confidence_Ca_o=scipy.stats.scoreatpercentile(all_output[:,9,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Ca_p=scipy.stats.scoreatpercentile(all_output[:,10,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_CO3_o=scipy.stats.scoreatpercentile(all_output[:,11,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_CO3_p=scipy.stats.scoreatpercentile(all_output[:,12,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_HCO3_o=scipy.stats.scoreatpercentile(all_output[:,13,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_HCO3_p=scipy.stats.scoreatpercentile(all_output[:,14,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)        
         
    pylab.subplot(3, 4, 3)
    pylab.plot(all_output[0,4,:],confidence_Ca_o[1],'r',label='Ca ocean')
    pylab.plot(all_output[0,4,:],confidence_Ca_p[1],'b',label='Ca pore')
    pylab.fill_between(all_output[0,4,:], confidence_Ca_o[0], confidence_Ca_o[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_Ca_p[0], confidence_Ca_p[2], color='red', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_HCO3_o[1],'k',label='HCO3 ocean')
    pylab.plot(all_output[0,4,:],confidence_HCO3_p[1],'g',label='HCO3 pore')
    pylab.fill_between(all_output[0,4,:],confidence_HCO3_o[0],confidence_HCO3_o[2], color='grey', alpha='0.4')
    pylab.fill_between(all_output[0,4,:],confidence_HCO3_p[0],confidence_HCO3_p[2], color='green', alpha='0.4')
    pylab.xlabel('Molality (mol/kg)')
    pylab.xlabel('Time (yr)')
    pylab.legend(loc=2)
    
    pylab.subplot(3, 4, 4)
    pylab.plot(all_output[0,4,:],confidence_CO3_o[1],'r',label='CO3 ocean')
    pylab.plot(all_output[0,4,:],confidence_CO3_p[1],'b',label='CO3 pore')
    pylab.fill_between(all_output[0,4,:], confidence_CO3_o[0], confidence_CO3_o[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_CO3_p[0], confidence_CO3_p[2], color='blue', alpha='0.4')
    pylab.xlabel('Time (yr)')
    pylab.ylabel('Molality (mol/kg)')
    pylab.legend()    
        
    pylab.subplot(3, 4, 5)
    confidence_omega_o=scipy.stats.scoreatpercentile(all_output[:,15,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_omega_p=scipy.stats.scoreatpercentile(all_output[:,16,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.plot(all_output[0,4,:],confidence_omega_o[1],'r',label='ocean')
    pylab.plot(all_output[0,4,:],confidence_omega_p[1],'b',label='pore space')
    pylab.fill_between(all_output[0,4,:], confidence_omega_o[0],confidence_omega_o[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_omega_p[0], confidence_omega_p[2], color='blue', alpha='0.4')
    pylab.legend(loc=2)
    pylab.ylabel('Saturation state')
    pylab.xlabel('Time (yr)')
    
    confidence_Tsurf=scipy.stats.scoreatpercentile(all_output[:,17,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Tdeep=scipy.stats.scoreatpercentile(all_output[:,18,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.subplot(3, 4, 6)
    pylab.fill_between(all_output[0,4,:], confidence_Tsurf[0], confidence_Tsurf[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_Tdeep[0], confidence_Tdeep[2], color='blue', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_Tsurf[1],'r',label='Surface')
    pylab.plot(all_output[0,4,:],confidence_Tdeep[1],'b',label='Deep')
    pylab.ylabel('Temperature (K)')
    pylab.xlabel('Time (yr)')
    pylab.legend(loc=2)
    
    confidence_Fd=scipy.stats.scoreatpercentile(all_output[:,19,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Fs=scipy.stats.scoreatpercentile(all_output[:,20,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_o=scipy.stats.scoreatpercentile(all_output[:,21,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_p=scipy.stats.scoreatpercentile(all_output[:,22,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.subplot(3, 4, 7)
    pylab.plot(all_output[0,4,:], confidence_Fs[1],'b',label='Continental weathering')
    pylab.plot(all_output[0,4,:],confidence_Prec_o[1],'g',label='ocean precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fs[0], confidence_Fs[2], color='blue', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_Prec_o[0], confidence_Prec_o[2], color='green', alpha='0.4')
    pylab.ylabel('Fluxes (mol C/yr)')
    pylab.xlabel('Time (yr)')     
    pylab.legend(loc=2)     
    
    pylab.subplot(3, 4, 8)
    pylab.plot(all_output[0,4,:],confidence_Fd[1],'r',label='Seafloor dissolution')
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1],'k',label='pore precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fd[0], confidence_Fd[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_Prec_p[0], confidence_Prec_p[2], color='grey', alpha='0.4')  

    spread_dist=[]
    for kk in range(0,10000):
        rate=1+numpy.random.uniform(0.2,1.5)
        beta_plot=numpy.random.uniform(0.0,1.0)
        precip_e=(numpy.random.uniform(-.75,0.75)+2.35)*10**12 
        spread_dist.append(rate**beta_plot*precip_e)   
    spread_dist=numpy.array(spread_dist)  
    [Slow,Smed,Shigh]=scipy.stats.scoreatpercentile(spread_dist,[5,50,95], interpolation_method='fraction',axis=0)    
    prec=Smed
    prec_er=numpy.array([[Smed-Slow],[Shigh-Smed]])  

    pylab.legend(loc=2)
    pylab.ylabel('Fluxes (mol C/yr)')
    pylab.xlabel('Time (yr)')
    
    pylab.subplot(3, 4, 10)
    confidence_DICo=scipy.stats.scoreatpercentile(all_output[:,0,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_ALKo=scipy.stats.scoreatpercentile(all_output[:,1,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_DICp=scipy.stats.scoreatpercentile(all_output[:,2,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_ALKp=scipy.stats.scoreatpercentile(all_output[:,3,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)

    pylab.plot(all_output[0,4,:],confidence_DICo[1],'r',label='DIC_o')
    pylab.plot(all_output[0,4,:],confidence_ALKo[1],'b',label='ALK_o')
    pylab.plot(all_output[0,4,:],confidence_DICp[1],'g',label='DIC_p')
    pylab.plot(all_output[0,4,:],confidence_ALKp[1],'k',label='ALK_p')
    
    pylab.fill_between(all_output[0,4,:], confidence_DICo[0], confidence_DICo[2], color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_ALKo[0], confidence_ALKo[2], color='blue', alpha='0.4')
    pylab.fill_between(all_output[0,4,:], confidence_DICp[0], confidence_DICp[2], color='green', alpha='0.4')
    
    pylab.xlabel('Time (yr)')
    pylab.ylabel('Molality (mol/kg)')
    pylab.legend()
    
    pylab.subplot(3, 4, 11)
    pylab.xlabel('relative change in silicate weathering')
    pylab.ylabel('Number')
    pylab.tight_layout()

    confidence_Fd=scipy.stats.scoreatpercentile(all_output[:,19,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Fs=scipy.stats.scoreatpercentile(all_output[:,20,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_o=scipy.stats.scoreatpercentile(all_output[:,21,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_p=scipy.stats.scoreatpercentile(all_output[:,22,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Volc=scipy.stats.scoreatpercentile(all_output[:,23,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_meta=scipy.stats.scoreatpercentile(all_output[:,24,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_sub=scipy.stats.scoreatpercentile(all_output[:,25,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_mantle_out=scipy.stats.scoreatpercentile(all_output[:,26,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)    

    confidence_Rcrust=scipy.stats.scoreatpercentile(all_output[:,27,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Rmantle=scipy.stats.scoreatpercentile(all_output[:,28,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Rorg=scipy.stats.scoreatpercentile(all_output[:,29,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_meta_carb=scipy.stats.scoreatpercentile(all_output[:,30,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_sub_carb=scipy.stats.scoreatpercentile(all_output[:,31,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_meta_org=scipy.stats.scoreatpercentile(all_output[:,32,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_sub_org=scipy.stats.scoreatpercentile(all_output[:,33,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_org_weath=scipy.stats.scoreatpercentile(all_output[:,34,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_org_burial=scipy.stats.scoreatpercentile(all_output[:,35,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_carb_weath=scipy.stats.scoreatpercentile(all_output[:,43,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_org_buried=scipy.stats.scoreatpercentile(all_output[:,44,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_carb_precip=scipy.stats.scoreatpercentile(all_output[:,45,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    
    forg_dist=all_output[:,35,:]/(all_output[:,35,:]+all_output[:,21,:]+all_output[:,22,:])
    confidence_forg_dist=scipy.stats.scoreatpercentile(forg_dist,[2.5,50,97.5], interpolation_method='fraction',axis=0)
    forg_apparent_dist=(all_output[:,38,:]-all_output[:,37,:])/(all_output[:,36,:]-all_output[:,37,:])
    confidence_forg_apparent=scipy.stats.scoreatpercentile(forg_apparent_dist,[2.5,50,97.5], interpolation_method='fraction',axis=0)

    forg_apparent_dist2=(all_output[:,38,:]-all_output[:,45,:])/(all_output[:,44,:]-all_output[:,45,:]) 
    confidence_forg_apparent2=scipy.stats.scoreatpercentile(forg_apparent_dist2,[2.5,50,97.5], interpolation_method='fraction',axis=0)

    forg_apparent_dist3=(-5.5 - all_output[:,45,:])/(all_output[:,44,:]-all_output[:,45,:])
    confidence_forg_apparent3=scipy.stats.scoreatpercentile(forg_apparent_dist3,[2.5,50,97.5], interpolation_method='fraction',axis=0)

    pylab.figure()
    pylab.subplot(2, 1, 1)
    pylab.plot(all_output[0,4,:],confidence_org_weath[1],'k',label='org weath',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_org_burial[1],'r',label='org burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1]+confidence_Prec_o[1],'c',label='carb burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1],'c--',label='pore burial')
    pylab.plot(all_output[0,4,:],confidence_Prec_o[1],'c*',label='ocean burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_org_weath[1]+confidence_carb_weath[1],'g*',label='Total weath',linewidth=2.5)
    pylab.ylabel('Flux (mol C/yr)')
    pylab.xlabel('Time (Ga)')
    pylab.legend()
    
    pylab.subplot(2, 1, 2)
    pylab.plot(all_output[0,4,:],confidence_forg_dist[1],'r',label='forg')
    pylab.fill_between(all_output[0,4,:], confidence_forg_dist[0], confidence_forg_dist[2], color='red', alpha='0.4',label='forg dist')
    pylab.plot(all_output[0,4,:],confidence_forg_apparent[1],'k--',label='Apparent forg - model')
    pylab.fill_between(all_output[0,4,:], confidence_forg_apparent[0], confidence_forg_apparent[2], color='grey', alpha='0.4',label='forg dist')
    LOWESS=numpy.loadtxt('forg_LOWESS.txt')
    pylab.plot(LOWESS[:,0]*-1e6,LOWESS[:,2],'b',label='Apparent forg - data')
    pylab.fill_between(LOWESS[:,0]*-1e6, LOWESS[:,1], LOWESS[:,3], color='blue', alpha='0.4',label='forg dist')
    pylab.ylabel('forg')
    pylab.xlabel('Time (Ga)')

    pylab.figure()
    pylab.subplot(2, 1, 1)
    pylab.plot(all_output[0,4,:],confidence_org_weath[1],'k',label='org weath',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_org_burial[1],'r',label='org burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1]+confidence_Prec_o[1],'c',label='carb burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1],'c--',label='pore burial')
    pylab.plot(all_output[0,4,:],confidence_Prec_o[1],'c*',label='ocean burial',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_org_weath[1]+confidence_carb_weath[1],'g*',label='Total weath',linewidth=2.5)
    pylab.ylabel('Flux (mol C/yr)')
    pylab.xlabel('Time (Ga)')
    pylab.legend()

    pylab.subplot(2, 1, 2)
    pylab.plot(all_output[0,4,:],confidence_forg_dist[1],'r',label='forg')
    pylab.fill_between(all_output[0,4,:], confidence_forg_dist[0], confidence_forg_dist[2], color='red', alpha='0.4',label='forg dist')
    pylab.plot(all_output[0,4,:],confidence_forg_apparent2[1],'k--',label='Apparent forg - model')
    pylab.fill_between(all_output[0,4,:], confidence_forg_apparent2[0], confidence_forg_apparent2[2], color='grey', alpha='0.4',label='forg dist')
    LOWESS=numpy.loadtxt('forg_LOWESS.txt')
    pylab.plot(LOWESS[:,0]*-1e6,LOWESS[:,2],'b',label='Apparent forg - data')
    pylab.fill_between(LOWESS[:,0]*-1e6, LOWESS[:,1], LOWESS[:,3], color='blue', alpha='0.4',label='forg dist')
    pylab.ylabel('forg')
    pylab.xlabel('Time (Ga)')

    carb_filt=numpy.ones(shape=(len(numpy.load("carb_200my_t.npy")),3))
    carb_filt[:,0]=numpy.load("carb_200my_t.npy") #time values for carbonates
    carb_filt[:,1]=numpy.load("carb_200my_v.npy") #isotope values for carbonates
    carb_filt[:,2]=numpy.load("carb_200my_var.npy") #isotope values for carbonates

    org_filt=numpy.ones(shape=(len(numpy.load("org_200my_t.npy")),3))
    org_filt[:,0]=numpy.load("org_200my_t.npy") #time values for organics
    org_filt[:,1]=numpy.load("org_200my_v.npy") #isotope values for organics
    org_filt[:,2]=numpy.load("org_200my_var.npy") #isotope values for organics

    carb_filt=numpy.ones(shape=(len(numpy.load("carb_RNG200my_t.npy")),3))
    carb_filt[:,0]=numpy.load("carb_RNG200my_t.npy") #time values for carbonates
    carb_filt[:,1]=numpy.load("carb_RNG200my_v.npy") #isotope values for carbonates
    carb_filt[:,2]=numpy.load("carb_RNG200my_var.npy") #isotope values for carbonates

    org_filt=numpy.ones(shape=(len(numpy.load("org_RNG200my_t.npy")),3))
    org_filt[:,0]=numpy.load("org_RNG200my_t.npy") #time values for organics
    org_filt[:,1]=numpy.load("org_RNG200my_v.npy") #isotope values for organics
    org_filt[:,2]=numpy.load("org_RNG200my_var.npy") #isotope values for organics

     #isotopes
    confidence_del_org=scipy.stats.scoreatpercentile(all_output[:,36,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_carb=scipy.stats.scoreatpercentile(all_output[:,37,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_mantle=scipy.stats.scoreatpercentile(all_output[:,38,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_AO=scipy.stats.scoreatpercentile(all_output[:,39,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_crust=scipy.stats.scoreatpercentile(all_output[:,40,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_out=scipy.stats.scoreatpercentile(all_output[:,41,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_del_input=scipy.stats.scoreatpercentile(all_output[:,42,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)     

    pylab.figure()
    pylab.plot(all_output[0,4,:],confidence_del_input[1],'k',label='input isotopes',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:], confidence_del_input[0], confidence_del_input[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_del_crust[1],'r',label='del crust',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_del_mantle[1], 'm',label='del mantle',linewidth=2.5)
    pylab.plot(all_output[0,4,:],confidence_del_out[1],'c',label='del output',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:], confidence_del_mantle[0], confidence_del_mantle[2], color='magenta', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_del_AO[1],'g',label='del AO',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:], confidence_del_AO[0], confidence_del_AO[2], color='green', alpha='0.4')

    pylab.figure()
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_input[1],'k',label='input isotopes',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_del_input[0], confidence_del_input[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_crust[1],'r',label='crustal isotopes',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_del_crust[0], confidence_del_crust[2], color='red', alpha='0.4')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_mantle[1], 'm',label='mantle isotopes',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_del_mantle[0], confidence_del_mantle[2], color='magenta', alpha='0.4')
    pylab.ylabel('Carbon isotope ratio')
    pylab.xlabel('Time (Ga)')
    pylab.legend(frameon=False)
    
    numpy.save("outputs_for_later_plotting3",all_output) 

    pylab.figure()
    org_crust_frac = all_output[:,29,:] /(all_output[:,29,:] + all_output[:,27,:])
    org_frac_conf = scipy.stats.scoreatpercentile(org_crust_frac,[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.plot(all_output[0,4,:]/-1e9,org_frac_conf[1],'k',label='organic crustal fraction',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, org_frac_conf[0], org_frac_conf[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_forg_dist[1],'r',label='True $f_{org}$ - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_forg_dist[0], confidence_forg_dist[2], color='red', alpha='0.4')

    pylab.figure()
    pylab.plot(carb_filt[:,0]*-1e6,carb_filt[:,1],'gx')
    pylab.plot(org_filt[:,0]*-1e6,org_filt[:,1],'ks')
    pylab.errorbar(org_filt[:,0]*-1e6,org_filt[:,1],yerr=org_filt[:,2],color='k',marker='o',linestyle="None")
    pylab.errorbar(carb_filt[:,0]*-1e6,carb_filt[:,1],yerr=carb_filt[:,2],color='g',marker='o',linestyle="None")

    pylab.plot(all_output[0,4,:],confidence_del_org[1],'k',label='org crust')
    pylab.fill_between(all_output[0,4,:], confidence_del_org[0], confidence_del_org[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_del_carb[1],'g',label='carb crust')
    pylab.fill_between(all_output[0,4,:], confidence_del_carb[0], confidence_del_carb[2], color='green', alpha='0.4')
 
    pylab.figure()
    pylab.plot(carb_filt[:,0]*-1e6,carb_filt[:,1],'gx')
    pylab.plot(org_filt[:,0]*-1e6,org_filt[:,1],'ks')
    pylab.errorbar(org_filt[:,0]*-1e6,org_filt[:,1],yerr=org_filt[:,2],color='k',marker='o',linestyle="None")
    pylab.errorbar(carb_filt[:,0]*-1e6,carb_filt[:,1],yerr=carb_filt[:,2],color='g',marker='o',linestyle="None")
 
    pylab.plot(all_output[0,4,:],confidence_del_org_buried[1],'k',label='org crust')
    pylab.fill_between(all_output[0,4,:], confidence_del_org_buried[0], confidence_del_org_buried[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:],confidence_del_carb_precip[1],'g',label='carb crust')
    pylab.fill_between(all_output[0,4,:], confidence_del_carb_precip[0], confidence_del_carb_precip[2], color='green', alpha='0.4')
    pylab.ylabel('Carbon isotope ratio')
    pylab.xlabel('Time (yrs)')


    pylab.figure()
    pylab.subplot(2, 1, 1)
    pylab.plot(all_output[0,4,:]/-1e9,confidence_forg_dist[1],'r',label='True $f_{org}$ - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_forg_dist[0], confidence_forg_dist[2], color='red', alpha='0.4')#,label='forg dist')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_forg_apparent2[1],'b--',label='Apparent $f_{org}$ - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_forg_apparent2[0], confidence_forg_apparent2[2], color='blue', alpha='0.4')#,label='forg dist')
    pylab.ylabel('Fractional organic burial, $f_{org}$')
    pylab.xlabel('Time (Ga)')
    pylab.legend(loc=3,frameon=False)

    pylab.subplot(2, 1, 2)
    pylab.plot(carb_filt[:,0]/1e3,carb_filt[:,1],'go',label='Carbonate burial - data')
    pylab.plot(org_filt[:,0]/1e3,org_filt[:,1],'ko',label='Organic burial - data')
    pylab.errorbar(carb_filt[:,0]/1e3,carb_filt[:,1],yerr=carb_filt[:,2],color='g',marker='o',linestyle="None")
    pylab.errorbar(org_filt[:,0]/1e3,org_filt[:,1],yerr=org_filt[:,2],color='k',marker='o',linestyle="None")
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_org_buried[1],'k',label='Organic burial - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_del_org_buried[0], confidence_del_org_buried[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_carb_precip[1],'g',label='Carbonate burial - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_del_carb_precip[0], confidence_del_carb_precip[2], color='green', alpha='0.4')
    pylab.xlabel('Time (Ga)')
    pylab.ylabel('Carbon isotopic ratio, $\delta^{13}$C')
    pylab.legend(loc=8,ncol=2,numpoints=1,frameon=False)
    pylab.ylim([-60,20])

    pylab.figure()
    pylab.subplot(2, 1, 1)
    pylab.plot(all_output[0,4,:]/-1e9,confidence_forg_dist[1],'r',label='True $f_{org}$ - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_forg_dist[0], confidence_forg_dist[2], color='red', alpha='0.4')#,label='forg dist')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_forg_apparent3[1],'b--',label='Apparent $f_{org}$ - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_forg_apparent3[0], confidence_forg_apparent3[2], color='blue', alpha='0.4')#,label='forg dist')
    pylab.ylabel('Fractional organic burial, $f_{org}$')
    pylab.xlabel('Time (Ga)')
    pylab.legend(loc=3,frameon=False)

    pylab.subplot(2, 1, 2)
    pylab.plot(carb_filt[:,0]/1e3,carb_filt[:,1],'go',label='Carbonate burial - data')
    pylab.plot(org_filt[:,0]/1e3,org_filt[:,1],'ko',label='Organic burial - data')
    pylab.errorbar(carb_filt[:,0]/1e3,carb_filt[:,1],yerr=carb_filt[:,2],color='g',marker='o',linestyle="None")
    pylab.errorbar(org_filt[:,0]/1e3,org_filt[:,1],yerr=org_filt[:,2],color='k',marker='o',linestyle="None")
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_org_buried[1],'k',label='Organic burial - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_del_org_buried[0], confidence_del_org_buried[2], color='grey', alpha='0.4')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_del_carb_precip[1],'g',label='Carbonate burial - model',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9,confidence_del_carb_precip[0], confidence_del_carb_precip[2], color='green', alpha='0.4')
    pylab.xlabel('Time (Ga)')
    pylab.ylabel('Carbon isotopic ratio, $\delta^{13}$C')
    pylab.legend(loc=8,ncol=2,numpoints=1,frameon=False)
    pylab.ylim([-60,20])


    ########################################
    ## FIGURE FOR PAPER
    ########################################

    all_output[0,4,:]=all_output[0,4,:]
    
    strt_lim=-4.5e9
    fin_lim=0.0
    pylab.figure(figsize=(30,15))
    pylab.subplot(3, 3, 1)
    pylab.plot(all_output[0,4,:]/-1e9,confidence_pH_o[1],'k',label='ocean')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_pH_o[0], confidence_pH_o[2], color='grey', alpha='0.4')
    matplotlib.pyplot.errorbar(all_output[0,4,999],true_mod_pH,yerr=er_true_mod_pH,color='k',marker='o')
    pylab.xlabel('Time (Ga)')
    pylab.ylabel('Ocean pH')

    pylab.subplot(3, 3, 2)
    pylab.semilogy(all_output[0,4,:]/-1e9,confidence_CO2o[1]/ppCO2,'k',label='Atmospheric CO2')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_CO2o[0]/ppCO2, confidence_CO2o[2]/ppCO2, color='grey', alpha='0.4')
    true_mod_CO2 = cons_dict['true_mod_CO2']/ppCO2
    er_CO2 = cons_dict['er_CO2']
    positive_er = ( -10**( numpy.log10(cons_dict['true_mod_CO2']) - er_CO2 ) + cons_dict['true_mod_CO2'] ) /ppCO2
    negative_er = ( 10**( numpy.log10(cons_dict['true_mod_CO2']) + er_CO2 ) - cons_dict['true_mod_CO2'] ) /ppCO2
    matplotlib.pyplot.errorbar(all_output[0,4,999],true_mod_CO2,yerr=numpy.transpose([[negative_er,positive_er]]),color='k',marker='o')
    pylab.xlabel('Time (Ga)')
    pylab.ylabel('Atmospheric pCO2 (ppm)')
   
    pylab.subplot(3, 3, 4)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Tsurf[0], confidence_Tsurf[2], color='grey', alpha='0.4')  
    pylab.plot(all_output[0,4,:]/-1e9,confidence_Tsurf[1],'k',label='Surface')
    true_mod_Temp = cons_dict['true_mod_Temp']
    er_modeTemp = cons_dict['er_modeTemp']
    pylab.errorbar(all_output[0,4,999],true_mod_Temp,yerr=er_modeTemp,color='k',marker='o',linestyle="None")

    pylab.ylabel('Mean Surface Temperature (K)')
    pylab.xlabel('Time (Ga)')

    pylab.subplot(3,3,5)
    pylab.plot(all_output[0,4,:]/-1e9, confidence_Fs[1]/1e12,'r',label='Cont. weathering')
    pylab.plot(all_output[0,4,:]/-1e9,confidence_Fd[1]/1e12,'k',label='Seafloor dissolution')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Fs[0]/1e12, confidence_Fs[2]/1e12, color='red', alpha='0.4')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Fd[0]/1e12, confidence_Fd[2]/1e12, color='grey', alpha='0.4')
    pylab.ylabel('Fluxes (Tmol C/yr)')
    pylab.xlabel('Time (Ga)')
    pylab.legend(loc=2)  
    
    pylab.subplot(3, 3, 6)
    pylab.plot(all_output[0,4,:]/-1e9,confidence_org_burial[1]/1e12,'k',label='org burial',linewidth=2.5)
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_org_burial[0]/1e12, confidence_org_burial[2]/1e12, color='grey', alpha='0.4')
    pylab.ylabel('Organic burial (Tmol C/yr)')
    pylab.xlabel('Time (Ga)')
    
    pylab.subplot(3, 3, 7)
    pO2_arr=all_output[:,46,:]
    pO2_confidence=scipy.stats.scoreatpercentile(pO2_arr,[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.semilogy(all_output[0,4,:]/-1e9, pO2_confidence[1],color='k')
    pylab.fill_between(all_output[0,4,:]/-1e9, pO2_confidence[0], pO2_confidence[2], color='grey', alpha='0.4')
    O2_source=all_output[:,35,:]+5.15e12*numpy.transpose((numpy.transpose(forg_dist)/forg_dist[:,999])) #5.15 is other O2 sources e.g.sulfide and iron oxide
    O2_sink=numpy.transpose(2.4e12*numpy.transpose(all_output[:,23,:])/numpy.transpose(all_output[:,23,999])) #+all_output[34,:,:]
    Kox=O2_source/O2_sink
    Kox_confidence=scipy.stats.scoreatpercentile(Kox,[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.ylabel('pO2 (PAL)')
    pylab.xlabel('Time (Ga)')
    
    pylab.subplot(3, 3, 3)
    pylab.plot(all_output[0,4,:]/-1e9,Kox_confidence[1],'k')
    pylab.fill_between(all_output[0,4,:]/-1e9, Kox_confidence[0], Kox_confidence[2],color='grey', alpha='0.4')
    pylab.ylabel('Relative oxygenation, K$_{oxy}$')
    pylab.xlabel('Time (Ga)')

    pylab.subplot(3, 3, 9)

    Avice_outgassing = cons_dict['Avice_outgassing']*confidence_mantle_out[1][999]/1e12 
    error_Avice_outgassing = cons_dict['error_Avice_outgassing']*confidence_mantle_out[1][999]/1e12  
    pylab.errorbar(all_output[0,4,194]/-1e9,Avice_outgassing,yerr=error_Avice_outgassing,color='r',marker='o',linestyle="None")
    pylab.plot(all_output[0,4,:]/-1e9,confidence_Volc[1]/1e12,'k',label='Total outgas.')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Volc[0]/1e12, confidence_Volc[2]/1e12, color='grey', alpha='0.4') 
    pylab.plot(all_output[0,4,:]/-1e9,confidence_mantle_out[1]/1e12,'r',label='Mantle outgas.')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_mantle_out[0]/1e12, confidence_mantle_out[2]/1e12, color='red', alpha='0.4') 
    pylab.ylabel('Outgassing (Tmol/yr)')
    pylab.xlabel('Time (Ga)')
    pylab.legend()

    pylab.subplot(3, 3, 8)
    pylab.semilogy(all_output[0,4,:]/-1e9,confidence_Rmantle[1],'k',label='Mantle')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Rmantle[0], confidence_Rmantle[2], color='grey', alpha='0.4') 
    pylab.semilogy(all_output[0,4,:]/-1e9,confidence_Rcrust[1],'r',label='Crustal carbonates')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Rcrust[0], confidence_Rcrust[2], color='red', alpha='0.4') 
    pylab.semilogy(all_output[0,4,:]/-1e9,confidence_Rorg[1],'g',label='Crustal organics')
    pylab.fill_between(all_output[0,4,:]/-1e9, confidence_Rorg[0], confidence_Rorg[2], color='green', alpha='0.4') 
    mantle_mod_true= cons_dict['mantle_mod_true'] 
    er_mantle_true = cons_dict['er_mantle_true'] 
    Rcrust_mod_true= cons_dict['Rcrust_mod_true'] 
    er_Rcrust_true = cons_dict['er_Rcrust_true'] 
    Rorg_mod_true= cons_dict['Rorg_mod_true'] 
    er_Rorg_true = cons_dict['er_Rorg_true'] 
    pylab.errorbar(all_output[0,4,999],mantle_mod_true,yerr=er_mantle_true,color='k',marker='o',linestyle="None")
    pylab.errorbar(all_output[0,4,999],Rcrust_mod_true,yerr=er_Rcrust_true,color='r',marker='o',linestyle="None")
    pylab.errorbar(all_output[0,4,999],Rorg_mod_true,yerr=er_Rorg_true,color='g',marker='o',linestyle="None")
    pylab.xlabel('Time (Ga)')
    pylab.ylabel('Reservoir size (mol C)')
    pylab.legend()
    pylab.show()
    
    pylab.tight_layout()
