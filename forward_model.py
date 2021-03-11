#! /usr/local/bin/python
import numpy
import pylab
import scipy.stats
from thermodynamic_variables import Sol_prod
from clim_model import clim_fun #cliamte model

def Precambrian(Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,del_mantle0,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,cons_dict,multplic,Mantle_temp_mod):


    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0,Ca_p,Ca_o,Rcrust0,Rorg0,k_meta0_carb,k_meta0_org,k_sub0_carb,k_sub0_org,F_orgB0,imply_forg0,org_weath0,carb_minus_org,carb_minus_AO
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2,ALK_o,ALK_p,k_meta0,k_sub0,k_mantle0,k_sub_arc0
    global mantle_mod,Rcrust_mod,Rorg_mod,F_meta_mod_carb,F_meta_mod_org
    global F_sub_mod_carb,F_sub_mod_org,ppCO2_mod,deep_grad

    T=18+273
    pK1=17.788 - .073104 *T - .0051087*35 + 1.1463*10**-4*T**2
    pK2=20.919 - .064209 *T - .011887*35 + 8.7313*10**-5*T**2
    H_CO2=pylab.exp(9345.17/T - 167.8108 + 23.3585 * pylab.log(T)+(.023517 - 2.3656*10**-4*T+4.7036*10**-7*T**2)*35)
    
    deep_grad = 1.0
    mantle_mod=cons_dict['mantle_mod_true']
    Rcrust_mod=cons_dict['Rcrust_mod_true'] 
    Rorg_mod=cons_dict['Rorg_mod_true'] 

    F_meta_mod_carb=F_meta_mod*Rcrust_mod/(Rcrust_mod+Rorg_mod)
    F_meta_mod_org=F_meta_mod*Rorg_mod/(Rcrust_mod+Rorg_mod)

    F_sub_mod_carb=F_sub_mod*Rcrust_mod/(Rcrust_mod+Rorg_mod)
    F_sub_mod_org=F_sub_mod*Rorg_mod/(Rcrust_mod+Rorg_mod)
    ppCO2_mod=0.000280
  
    ############################################################################
    ############################################################################
    ## Find model initial conditions 
    Total_modern_inputs = F_out_mantle_mod + F_meta_mod + F_sub_mod + F_carbw + org_weath_mod0 +thermo_weath
    Modern_carbonate_burial = Total_modern_inputs * (1 - j1*j2*j3)
    k_w = Modern_carbonate_burial - 0.5e12 - F_carbw 

    if k_w < 0:
        print(k_w/1e12)
        print ("NEW PROBLEM: modern silciate weathering implied negative")
        return  -1

    Q0=(1-abs(4.1e9)/(4.5e9))**-n_out

    Mantle_temp0= Mantle_temp_mod*(Q0)**0.1
    Ccrit0 = (3.125e-3*80.0+835.5 - 285.0)/(Mantle_temp0-285.0)
    carb_eff_modifier0 = numpy.min([1.0,(Ccrit0 - 0.3)/(0.6-0.3)])
    if carb_eff_modifier0 < 0.0:
        carb_eff_modifier0 = 0.0

    Rcrust0=R_carb_00 
    Rorg0=0.0
    Rmantle0=R_mantle_00
    
    k_meta0_carb=F_meta_mod_carb*(Q0**mm)*(Rcrust0/Rcrust_mod)
    k_meta0_org=F_meta_mod_org*(Q0**mm)*(Rorg0/Rorg_mod)
    
    
    k_sub0_carb=F_sub_mod_carb*(Q0**beta)*(Rcrust0/Rcrust_mod)   
    k_sub0_org=F_sub_mod_org*(Q0**beta)*(Rorg0/Rorg_mod)
    
    F_mantle_out0=F_out_mantle_mod*(Q0**mm)*(Rmantle0/mantle_mod)
    F_outgass0=F_mantle_out0+F_meta_mod_carb*(Q0**mm)*(Rcrust0/Rcrust_mod) + F_meta_mod_org*(Q0**mm)*(Rorg0/Rorg_mod)+((1-carb_eff_modifier0)*F_sub_mod_carb*(Q0**beta)*(Rcrust0/Rcrust_mod)+e1*F_sub_mod_org*(Q0**beta)*(Rorg0/Rorg_mod))
    
    lum_vary="y" 
    landf_vary="y"
    
    ### Initial conditions for modern atmosphere and ocean
    pH_o=pH_00 
    ppCO2_o=ppCO2_00 
    CO2_aq_o=H_CO2 * ppCO2_o
    hco3_o=(10**-pK1) * CO2_aq_o / (10**-pH_o)
    co3_o=(10**-pK2) * hco3_o / (10**-pH_o)
    DIC_o=co3_o+hco3_o+CO2_aq_o
    ALK_o=2*co3_o+hco3_o
    Ca_o=0.10100278
    salt = ALK_o-2*Ca_o 
    
    T_surface_mod=clim_fun(ppCO2_mod,1/(1+0.4*(abs(0)/4.6e9)))
    T_surface0=clim_fun(ppCO2_o,1/(1+0.4*(abs(4.1e9)/4.6e9)))
    T_surface_diff=T_surface0-T_surface_mod
    interc=274.037-deep_grad*285.44411576277344
    buffer_T=deep_grad*T_surface0+interc
    buffer_Ti=numpy.max([deep_grad*T_surface0+interc,271.15])
    buffer_T=numpy.min([buffer_Ti,T_surface0])
    Pore_mod=9.0            

    biology_mod=1-1/(1/(1-Bio_f)+numpy.exp(-5*(abs(0.0)-1e9)/1e9))
    bio_Archean = 1-1/(1/(1-Bio_f)+numpy.exp(-5*(abs(4.1e9)-1e9)/1e9))
    biology_modifier=(bio_Archean/bio_Archean)/(biology_mod/bio_Archean)
    if landf_vary=="y":
        land=lf(4.1e9,lfrac,growth_timing)#/lf(4.1e9,lfrac,growth_timing)
        land_mod=lf(0,lfrac,growth_timing)#/lf(4.1e9,lfrac,growth_timing)
    else:
        land=1
        land_mod=1
                                                                                 
    carb_weath0=(Rcrust0/Rcrust_mod)*F_carbw*biology_modifier*(ppCO2_00/ppCO2_mod)**(climp)*numpy.exp((T_surface_diff)/tdep_weath)*(land/land_mod)    
 
    #initialize isotopes
    delMantle=del_mantle0 
    org_weath0=0.0  
    delOrg=del_org0
    delCarb=del_carb0 
    del_outgas0=(delMantle*F_mantle_out0+(k_meta0_carb+k_sub0_carb*(1- carb_eff_modifier0))*delCarb+(k_meta0_org+k_sub0_org*e1)*delOrg)/(F_outgass0)
    del_input0=(del_outgas0*F_outgass0+delCarb*carb_weath0 +org_weath0*delOrg)/(F_outgass0+carb_weath0 +org_weath0)
    carb_minus_AO=0 
    delAO=delCarb-carb_minus_AO
    carb_minus_org=delCarb-delOrg
    imply_forg0=0.0
    F_orgB0=0.0   
    Mo=1.35e21
    buffer_T_mod=deep_grad*T_surface_mod+interc
    buffer_Ti_mod=numpy.max([deep_grad*T_surface_mod+interc,271.15])
    buffer_T_mod=numpy.min([buffer_Ti_mod,T_surface_mod])
    F_sil0=k_w*biology_modifier*(land/land_mod)*(ppCO2_o/ppCO2_mod)**climp*numpy.exp((T_surface_diff)/tdep_weath)
    F_diss0=(F_outgass0+org_weath0)-F_sil0-F_orgB0
    if F_diss0<0:
        print(F_outgass0/1e12,F_diss0/1e12,F_sil0/1e12)
        print ("PROBLEM: negative silciate weathering because outgassing too small")
        return  -1
    k_r = F_diss0 / ((Q0**beta)*2.88*10**-14*10**(-coef_for_diss*pH_o)*numpy.exp(-Ebas/(8.314*(buffer_T+Pore_mod))))   
    Fprecip_p=F_diss0
                         
    # Solving system of equations for pore space properties
    DIC_p=DIC_o
    ALK_p=ALK_o
    Fprecip_o=(F_outgass0+carb_weath0+org_weath0-F_orgB0)-Fprecip_p    
    if DIC_p==0:
        print ("WARNING: NO PHYSIAL STEADY STATE EXISTS! DIC is negative")   
    H_p=10**(-pH_00)
    pH_p=-pylab.log10(H_p)
    CO2aq_p=CO2_aq_o
    co3_p=co3_o 
    hco3_p=hco3_o
    Ca_p = Ca_o
    omega_p=Ca_p *co3_p/Sol_prod(buffer_T+Pore_mod)
    omega_o=Ca_o *co3_o/Sol_prod(T_surface0)
    ## All initital conditions have been calculated
    ############################################################################
    ############################################################################
    
    ############ Calculate proportionality constants ###########################
    # Given modern precipitation fluxes and saturation states, calculate proportionality constants
    k_c=Fprecip_p/(omega_p )**n
    k_o=Fprecip_o/(omega_o )**n

    ## Partition ocean carbonate sink into shelf precipitation and carbonate precipitation
    frac_pel=0
    if lfrac<1.0:
        ko1=(1-frac_pel)*Fprecip_o/((omega_o)**n)
    else:
        ko1=(1-frac_pel)*Fprecip_o/((land/land_mod)*(omega_o)**n)
    ko2=frac_pel*Fprecip_o/(omega_o**2.84)
    
    ## Define time interval and run model
    time_ar=numpy.linspace(-4.1e9,0,1000) #negative works with negative DEs
    import time
    ta = time.time()
    [out,mes]=scipy.integrate.odeint(system_of_equations, [DIC_o+1.8e20/Mo*ppCO2_o,ALK_o,DIC_p,ALK_p,Rcrust0,Rmantle0,1e14,delOrg,delCarb,delMantle,delAO], time_ar, args=(Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,multplic,Mantle_temp_mod),full_output=1)#,mxstep=100000)
    tb = time.time()
    print ("b-a",tb-ta)
       
    ## Create output vectors for variables of interest
    pH_array_o=0*time_ar
    CO2_array_o=0*time_ar
    pH_array_p=0*time_ar
    CO2_array_p=0*time_ar
    Ca_array_o=0*time_ar
    Ca_array_p=0*time_ar
    CO3_array_o=0*time_ar
    CO3_array_p=0*time_ar
    HCO3_array_o=0*time_ar
    HCO3_array_p=0*time_ar
    omega_o=0*time_ar
    omega_p=0*time_ar
    Tsurf_array=0*time_ar
    Tdeep_array=0*time_ar
    Fd_array=0*time_ar
    Fs_array=0*time_ar
    Precip_ocean_array=0*time_ar
    Precip_pore_array=0*time_ar
    volc_array=0*time_ar
    carbw_ar=0*time_ar
    co2_aq_o=0*time_ar
    co2_aq_p=0*time_ar
    F_meta_ar=0*time_ar
    F_sub_ar=0*time_ar
    F_out_mantle_ar=0*time_ar
    F_meta_carb_ar=0*time_ar
    F_sub_carb_ar=0*time_ar
    F_meta_org_ar=0*time_ar
    F_sub_org_ar=0*time_ar
    org_weath_ar=0*time_ar
    F_orgB_ar=0*time_ar
    del_crust_ar=0*time_ar
    del_out_ar=0*time_ar
    del_inputs_ar=0*time_ar
    del_org_buried_ar=0*time_ar
    del_carb_precip_ar=0*time_ar
    pO2_ar=0*time_ar
    reduced_early_ar=0*time_ar
    F_sub_arc_ar = 0*time_ar
    del_arc_volc_ar = 0*time_ar
    del_meta_volc_ar = 0*time_ar
    del_sub_ar = 0*time_ar
    one_minus_carb_eff_modifier_ar = 0*time_ar
    Kox_array = 0*time_ar
    
    for i in range(0,len(time_ar)):  ## Fill output vectors
        pl1=carbon_cycle([out[i,0],out[i,1],out[i,2],out[i,3],out[i,4],out[i,5],out[i,6],out[i,7],out[i,8],out[i,9],out[i,10]],time_ar[i],Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,multplic,Mantle_temp_mod) 
        pH_array_o[i]=pl1[0]
        co2_aq_o[i]=pl1[3]
        CO2_array_o[i]=pl1[4]
        pH_array_p[i]=pl1[19]
        CO2_array_p[i]=pl1[20]
        Ca_array_o[i]=pl1[5]
        Ca_array_p[i]=pl1[21]
        CO3_array_o[i]=pl1[1]
        CO3_array_p[i]=pl1[16]
        HCO3_array_o[i]=pl1[2]
        HCO3_array_p[i]=pl1[17]
        co2_aq_p[i]=pl1[18]
        omega_o[i]=pl1[6]
        omega_p[i]=pl1[13]
        Tsurf_array[i]=pl1[7] 
        Tdeep_array[i]=pl1[8]
        Fd_array[i]=pl1[9]  
        Fs_array[i]=pl1[10]
        Precip_ocean_array[i]=pl1[11]
        Precip_pore_array[i]=pl1[12]
        carbw_ar[i]=pl1[15]
        volc_array[i]=pl1[14]
        F_meta_ar[i]=pl1[26]
        F_sub_ar[i]=pl1[27]
        F_out_mantle_ar[i]=pl1[28] 
        F_meta_carb_ar[i]=pl1[29]
        F_sub_carb_ar[i]=pl1[30]
        F_meta_org_ar[i]=pl1[31]
        F_sub_org_ar[i]=pl1[32]
        org_weath_ar[i]=pl1[33]
        F_orgB_ar[i]=pl1[34]
        del_crust_ar[i]=pl1[37]
        del_out_ar[i]=pl1[38]        
        del_inputs_ar[i]=pl1[39]
        del_org_buried_ar[i]=pl1[35]
        del_carb_precip_ar[i]=pl1[36]
        pO2_ar[i]=pl1[40]
        reduced_early_ar[i] = pl1[41]
        F_sub_arc_ar[i] = pl1[42]
        del_arc_volc_ar[i] = pl1[43]
        del_meta_volc_ar[i] = pl1[44]
        del_sub_ar[i] = pl1[45]
        one_minus_carb_eff_modifier_ar[i] = pl1[46]
        Kox_array[i] = pl1[47]

    imbalance=0
    

    return [out[:,0],out[:,1],out[:,2],out[:,3],time_ar,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array,volc_array,F_meta_ar,F_sub_ar,F_out_mantle_ar,out[:,4],out[:,5],out[:,6],F_meta_carb_ar,F_sub_carb_ar,F_meta_org_ar,F_sub_org_ar,org_weath_ar,F_orgB_ar,out[:,7],out[:,8],out[:,9],out[:,10],del_crust_ar,del_out_ar,del_inputs_ar,carbw_ar,del_org_buried_ar,del_carb_precip_ar,pO2_ar,reduced_early_ar, F_sub_arc_ar,del_arc_volc_ar, del_meta_volc_ar, del_sub_ar,  one_minus_carb_eff_modifier_ar,Kox_array],imbalance 
   
    
### System of equations for reservoir and isotope evolution:    
def system_of_equations(y,t0,Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,multplic,Mantle_temp_mod):
   
    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0,Rcrust0,Rorg0,k_meta0_carb,k_meta0_org,k_sub0_carb,k_sub0_org,F_orgB0,imply_forg0,org_weath0,carb_minus_org,carb_minus_AO
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2,k_meta0,k_sub0,k_mantle0,k_sub_arc0
    
    s=1.8e20/Mo #correction factor for mass balance 
    
    
    [EE_pH_o,EE_co3_o,EE_hco3_o, EE_co2aq_o,EE_ppCO2_o,EE_Ca_o,EE_omega_o,T_surface,buffer_T,F_diss,F_silic,Precip_ocean,Precip_pore,EE_omega_p,F_outg,carb_weath,EE_co3_p,EE_hco3_p, EE_co2aq_p,EE_pH_p,EE_ppCO2_p,EE_Ca_p,EE_H_o,EE_H_p,T_surface_diff,F_outg,F_meta,F_sub,F_out_mantle,F_meta_carb,F_sub_carb,F_meta_org,F_sub_org,org_weath,F_orgB,del_org_buried,del_carb_precip,del_crust,del_out,del_inputs,pO2,reduced_early, F_sub_arc,del_arc_volc,del_meta_volc, del_sub,one_minus_carb_eff_modifier,Kox]=carbon_cycle([y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10]], t0, Bio_f, ppCO2_00, n, climp, tdep_weath, F_meta_mod, F_out_mantle_mod, pH_00, lfrac, R_carb_00, del_carb0,F_carbw,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,multplic,Mantle_temp_mod)
    DICo=y[0]-EE_ppCO2_o*s 
    
    # Calculate time derivatives from current state (equation 6).
    dy0_dt=F_outg/Mo+carb_weath/Mo-Precip_ocean/Mo+org_weath/Mo-F_orgB/Mo- Precip_pore/Mo # time derivative of atmosphere-ocean carbon abundance

    dy1_dt=2*carb_weath/Mo+2*F_silic/Mo-2*Precip_ocean/Mo+2*F_diss/Mo-2*Precip_pore/Mo #Time derivative of ocean alkalinity

    dy2_dt=0.0 #Time derivative of pore space carbon abundance

    dy3_dt=0.0 #Time derivative of pore space alkalinity

    dy4_dt=Precip_ocean+Precip_pore-carb_weath-F_meta_carb-F_sub_carb #time rate change continental carbonates

    Mantle_temp= Mantle_temp_mod*((1-abs(t0)/(4.5e9))**-n_out)**0.1
    Ccrit = (3.125e-3*80.0+835.5 - 285.0)/(Mantle_temp-285.0)
    carb_eff_modifier = numpy.min([1.0,(Ccrit - 0.3)/(0.6-0.3)])
    if carb_eff_modifier < 0.0:
        carb_eff_modifier = 0.0

    dy5_dt=F_sub_org*(1-e1)+F_sub_carb*carb_eff_modifier-F_out_mantle #mantle reservoir evolution

    dy6_dt=F_orgB-org_weath-F_sub_org-F_meta_org # organic crustal reservoir
    
    dy7_dt=-dy6_dt*y[7]/y[6]+(1/y[6])*(F_orgB*del_org_buried-org_weath*y[7]-F_sub_org*y[7]-F_meta_org*y[7]) # change in d13C organic reseroir
    dy8_dt=-dy4_dt*y[8]/y[4]+(1/y[4])*(Precip_ocean*del_carb_precip+Precip_pore*del_carb_precip-carb_weath*y[8]-F_meta_carb*y[8]-F_sub_carb*y[8]) # change in d13C carb reseroir
    dy9_dt=-dy5_dt*y[9]/y[5]+(1/y[5])*(F_sub_org*y[7]*(1-e1)+F_sub_carb*y[8]*carb_eff_modifier-F_out_mantle*y[9]) # change in d13C mantle reseroir
    dyAO_dt=Mo*dy0_dt
    dy10_dt=-dyAO_dt*y[10]/(y[0]*Mo)+(1/(y[0]*Mo))*(F_outg*del_out+carb_weath*y[8]-Precip_ocean*del_carb_precip+org_weath*y[7]-F_orgB*del_org_buried-Precip_pore*del_carb_precip) #change in d13C AO reservoir
        
    return [dy0_dt,dy1_dt,dy2_dt,dy3_dt,dy4_dt,dy5_dt,dy6_dt,dy7_dt,dy8_dt,dy9_dt,dy10_dt]           

def lf (t0,lfrac,growth_timing): #Calculate continental land fraction
    land_frac=numpy.max([0.0,1-1/(lfrac+numpy.exp(-10*(abs(t0)/1e9-growth_timing)))])
    return land_frac

## Calculate carbon fluxes and ocean chemistry
def carbon_cycle (y,t0,Bio_f,ppCO2_00,n,climp,tdep_weath,F_meta_mod,F_out_mantle_mod,pH_00,lfrac,R_carb_00,del_carb0,F_carbw,R_mantle_00,F_sub_mod,coef_for_diss,beta,n_out,mm,growth_timing,new_add_Ca,Ebas,e1,protero_O2,del_org0,e2,j1,j2,j3,org_weath_mod0,thermo_weath,multplic,Mantle_temp_mod):
 
    Ca_cofactor=new_add_Ca/5000.0*((1-abs(t0)/(4.5e9))**-0.5-1) #Redundant 
    
    global pK1,pK2,H_CO2,salt,Pore_mod,T_surface0,ALK_o,ALK_p,Ca_p,Ca_o,Rcrust0,Rorg0,k_meta0_carb,k_meta0_org,k_sub0_carb,k_sub0_org,F_orgB0,imply_forg0,org_weath0,carb_minus_org,carb_minus_AO
    global landf_vary,lum_vary,salt0,k_w,k_c,k_o,k_r,Mp,ppCO2_o,Mo,ko1,ko2,k_meta0,k_sub0,k_mantle0,k_sub_arc0
    global mantle_mod,Rcrust_mod,Rorg_mod,F_meta_mod_carb,F_meta_mod_org
    global F_sub_mod_carb,F_sub_mod_org,ppCO2_mod,deep_grad

    s=1.8e20/Mo
    [c,d]=pylab.roots([y[1]/(10**-pK2*10**-pK1)*(1+s/H_CO2),(y[1]-y[0])/(10**-pK2),(y[1]-2*y[0])])
    EE_H_o=numpy.max([c ,d])
    EE_pH_o=-pylab.log10(EE_H_o)
    EE_co3_o=y[1]/(2+EE_H_o/(10**-pK2))
    EE_hco3_o=y[1]-2*EE_co3_o
    EE_co2aq_o=( EE_hco3_o*EE_H_o/(10**-pK1) )
    EE_ppCO2_o = EE_co2aq_o /H_CO2
    EE_Ca_o = Ca_o+0.5*(y[1]-ALK_o)+Ca_cofactor

    #pore space
    EE_H_p=EE_H_o
    EE_pH_p=EE_pH_o
    EE_co3_p=EE_co3_o
    EE_hco3_p=EE_hco3_o
    EE_co2aq_p= EE_co2aq_o
    EE_ppCO2_p = EE_ppCO2_o 
    EE_Ca_p = EE_Ca_o 

    bio_now = 1-1/(1/(1-Bio_f)+numpy.exp(-5*(abs(t0)-1e9)/1e9))
    bio_Archean = 1-1/(1/(1-Bio_f)+numpy.exp(-5*(abs(4.1e9)-1e9)/1e9))
    biology_modern =   1-1/(1/(1-Bio_f)+numpy.exp(-5*(abs(0.0)-1e9)/1e9))
    biology = bio_now/bio_Archean
    biology_mod = biology_modern/bio_Archean

    Q0=(1-abs(4.1e9)/(4.5e9))**-n_out
    Q=(1-abs(t0)/(4.5e9))**-n_out
    spread= Q**beta     

    Mantle_temp= Mantle_temp_mod*(Q)**0.1
    Ccrit = (3.125e-3*80.0+835.5 - 285.0)/(Mantle_temp-285.0)
    carb_eff_modifier = numpy.min([1.0,(Ccrit - 0.3)/(0.6-0.3)])
    if carb_eff_modifier < 0.0:
        carb_eff_modifier = 0.0
            
    F_out_mantle = F_out_mantle_mod*(Q**mm)*(y[5]/mantle_mod)
    F_meta_carb=F_meta_mod_carb*(Q**mm)*(y[4]/Rcrust_mod)
    F_meta_org=F_meta_mod_org*(Q**mm)*(y[6]/Rorg_mod)
    F_meta=F_meta_carb+F_meta_org
    F_sub_carb = F_sub_mod_carb*(Q**beta)*(y[4]/Rcrust_mod)
    F_sub_org = F_sub_mod_org*(Q**beta)*(y[6]/Rorg_mod)
    
    F_sub= F_sub_carb+F_sub_org
    F_sub_arc = F_sub_carb*(1- carb_eff_modifier)+F_sub_org*e1       
    F_outg=F_out_mantle+F_meta+F_sub_arc

    #isotopes
    del_sub = (F_sub_carb*y[8]+F_sub_org*y[7])/F_sub_arc
    del_arc_volc = (F_sub_carb*(1 - carb_eff_modifier)*y[8]+e1*F_sub_org*y[7])/F_sub_arc
    del_meta_volc = (F_meta_carb*y[8]+F_meta_org*y[7])/F_meta
    del_crust=(y[4]*y[8]+y[6]*y[7])/(y[4]+y[6])
    del_out=(F_out_mantle*y[9]+F_meta_carb*y[8]+F_meta_org*y[7]+F_sub_carb*(1 - carb_eff_modifier)*y[8]+(e1)*F_sub_org*y[7])/F_outg
    carb_minus_org=0
    carb_minus_AO=0
    if abs(t0)<4.0e9:
        carb_minus_org=28.0
        carb_minus_AO=2
    del_carb_precip=y[10]+carb_minus_AO 
    del_org_buried=del_carb_precip-carb_minus_org     
    
    if landf_vary=="y":
        land=lf(t0,lfrac,growth_timing)
        land_mod=lf(0,lfrac,growth_timing)
    else:
        land=1
        land_mod=1

    if lum_vary=="y":
        L_over_Lo = 1/(1+0.4*(abs(t0)/4.6e9)) 
        T_surface=clim_fun(EE_ppCO2_o,L_over_Lo)
    else:
        L_over_Lo = 1/(1+0.4*(abs(4.1e9)/4.6e9))
        T_surface=clim_fun(EE_ppCO2_o,L_over_Lo)

    T_surface_mod=clim_fun(ppCO2_mod,1/(1+0.4*(abs(0)/4.6e9)))
    T_surface_diff_mod=T_surface-T_surface_mod
    T_surface_diff=T_surface-T_surface0
    interc=274.037-deep_grad*285.44411576277344
    buffer_T=deep_grad*T_surface+interc
    buffer_Ti=numpy.max([deep_grad*T_surface+interc,271.15])
    buffer_T=numpy.min([buffer_Ti,T_surface])

    EE_omega_o=EE_Ca_o *EE_co3_o/(Sol_prod(T_surface))
    EE_omega_p=EE_Ca_p *EE_co3_p/(Sol_prod(buffer_T+Pore_mod))    
                   
    #ocean carbonate sink
    Carb_sink=ko1*(land/land_mod)*(EE_omega_o)**n+ko2*(EE_omega_o**2.84)
    Precip_ocean=Carb_sink
    Precip_pore=k_c*(EE_omega_p )**n   

    # Carbonate weathering flux 
    carb_weath=(y[4]/Rcrust_mod)*F_carbw*(biology/biology_mod)*(EE_ppCO2_o/ppCO2_mod)**(climp)*numpy.exp((T_surface_diff_mod)/tdep_weath)*(land/land_mod)
    F_diss = k_r *(2.88*10**-14*10**(-coef_for_diss*EE_pH_p)*numpy.exp(-Ebas/(8.314*(buffer_T+Pore_mod)))) * spread
    F_silic=k_w*(biology/biology_mod)*(land/land_mod)*(EE_ppCO2_o/ppCO2_mod)**climp*numpy.exp((T_surface_diff_mod)/tdep_weath)
   
    aa = -0.5*numpy.log10(protero_O2)
    bb= 4.5 + 0.5*numpy.log10(protero_O2)
    cc = -4.5
    pO2 = 10** ( aa * numpy.tanh(-20*(abs(t0)-0.55e9)/1e9)+bb * numpy.tanh(-20*(abs(t0)-2.45e9)/1e9) + cc)
    
    org_weath = org_weath_mod0*(y[6]/Rorg_mod)*(land/land_mod)*(pO2)**0.30809  +thermo_weath*(y[6]/Rorg_mod) ## This is running in retrieval

    total_carb_inputs = (carb_weath+org_weath+F_outg)
    if (abs(t0)>4.0e9):
        F_orgB=0*total_carb_inputs
    if (abs(t0)<4.0e9)and(abs(t0)>3.9e9):
        F_orgB=j1*total_carb_inputs*(-abs(t0)+4.0e9)/0.1e9
    if (abs(t0)<3.9e9)and(abs(t0)>2.5e9):
        F_orgB=j1*total_carb_inputs
    if (abs(t0)<2.5e9)and(abs(t0)>2.4e9):
        inbetween = (abs(t0)-2.4e9)/0.1e9 
        F_orgB=j1*total_carb_inputs*inbetween + (1-inbetween)*j2*j1*total_carb_inputs
    if (abs(t0)<2.4e9)and(abs(t0)>0.6e9):
        F_orgB=j2*j1*total_carb_inputs
    if (abs(t0)>0.5e9)and(abs(t0)<0.6e9):
        inbetween = (abs(t0)-0.5e9)/0.1e9 
        F_orgB=j1*j2*total_carb_inputs*inbetween + (1-inbetween)*j3*j2*j1*total_carb_inputs
    if (abs(t0)<0.5e9):
        F_orgB=j3*j2*j1*total_carb_inputs
    
    if (abs(t0)>4.0e9):
        tempforg=0
    if (abs(t0)<4.0e9)and(abs(t0)>3.9e9): 
        tempforg=j1*(-abs(t0)+4.0e9)/0.1e9
    if (abs(t0)<3.9e9)and(abs(t0)>2.5e9):
        tempforg=j1
    if (abs(t0)<2.5e9)and(abs(t0)>2.4e9):
        inbetween = (abs(t0)-2.4e9)/0.1e9 
        tempforg=j1*inbetween + (1-inbetween)*j2*j1
    if (abs(t0)<2.4e9)and(abs(t0)>0.6e9):
        tempforg=j2*j1
    if (abs(t0)>0.5e9)and(abs(t0)<0.6e9):
        inbetween = (abs(t0)-0.5e9)/0.1e9 
        tempforg=j1*j2*inbetween + (1-inbetween)*j3*j2*j1
    if (abs(t0)<0.5e9):
        tempforg=j3*j2*j1
  
    current_volc = F_out_mantle_mod + F_meta_mod + F_sub_mod
    multiplicative = 1 + multplic*abs(t0)/4.1e9
    O2_fast_sink = multiplicative*2.4e12* F_outg/current_volc 
    model_modern_forg = j1*j2*j3
    
    reduced_early= (F_orgB + 5.15e12 * tempforg/model_modern_forg - org_weath -  O2_fast_sink)       
    del_inputs=(F_outg*del_out+org_weath*y[7]+carb_weath*y[8])/(F_outg+org_weath+carb_weath)
    O2_source=F_orgB+5.15e12*(tempforg/model_modern_forg)
    Kox= O2_source/O2_fast_sink

    return [EE_pH_o,EE_co3_o,EE_hco3_o, EE_co2aq_o,EE_ppCO2_o,EE_Ca_o,EE_omega_o,T_surface,buffer_T,F_diss,F_silic,Precip_ocean,Precip_pore,EE_omega_p,F_outg,carb_weath,EE_co3_p,EE_hco3_p, EE_co2aq_p,EE_pH_p,EE_ppCO2_p,EE_Ca_p,EE_H_o,EE_H_p,T_surface_diff,F_outg,F_meta,F_sub,F_out_mantle,F_meta_carb,F_sub_carb,F_meta_org,F_sub_org,org_weath,F_orgB,del_org_buried,del_carb_precip,del_crust,del_out,del_inputs,pO2,reduced_early,F_sub_arc,del_arc_volc,del_meta_volc,del_sub,1- carb_eff_modifier,Kox] 
