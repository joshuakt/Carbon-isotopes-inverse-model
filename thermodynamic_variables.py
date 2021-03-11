import pylab

def Sol_prod(T): # from appendix G Pilson book, reproduces table G.1 for salinity=35
    bo=-.77712
    b1=0.0028426
    b2=178.34
    co=-.07711
    do=.0041249
    S=35
    (bo+b1*T+b2/T)*S**0.5+co*S+do*S**1.5
    logK0=-171.9065-0.077993*T+2839.319/T+71.595*pylab.log10(T)
    
    logK=logK0+(bo+b1*T+b2/T)*S**0.5+co*S+do*S**1.5
    return 10**logK


