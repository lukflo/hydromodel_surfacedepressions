import numpy as np
import vanGenuchten as vg
    
def subsurfN(pp,cnic,cnhic, cno2):

    mu = pp[0]
    sig = pp[1]
    
    mu2 = pp[5]
    sig2 = pp[6]
    
    mu3 = pp[8]
    sig3 = pp[9]
    
    a6 = pp[15]
    minerr = pp[7]

    mu_NH4 = pp[2] # 1/day
    mu_NO2 = pp[3]
    mu_NO3 = pp[4]
    
    rho_b = 1.3 # g/cm^3
    rho_b = rho_b*100*100*100/1000 # kg/m^3
    
    KOC_NH4 = 0.000034 
    KOC_NH4 = 0.0

    K_NH4 = 2.0 # grams per mass^3 in fluid phase == mg/l
    K_NO2 = 2.0
    K_NO3 = 2.0

    timeings = 1079
    
    delt = 1/96
    SimTime = 900*timeings/86400
    num_time = int(SimTime/delt) + 1
    
    pset = np.loadtxt('input_param.txt')
    
    dz = 0.10
    ProfileDepth = 3
    z = np.arange(dz/2.0, ProfileDepth, dz)
    spac = dz*np.ones(len(z)+0)

    PROF = np.load('Profiles.npy')
    VELS = np.load('Velocities.npy')
    SINK = np.load('Tiles.npy')
    
    pr = [PROF[0], PROF[1], PROF[2], PROF[3]]
    vels = [VELS[0], VELS[1], VELS[2], VELS[3]]
    sink_pr = [86400*SINK[0],86400*SINK[1],86400*SINK[2],86400*SINK[3]]


    lp = 30 # length of profile
    npr = len(pr) # number of profiles
    nt = len(pr[0]) # number of times
    
    MatricPot = np.empty((int(lp),int(npr),int(nt)))
    WaterCon = np.empty((int(lp),int(npr),int(nt)))
    Velocity = np.empty((int(lp+1),int(npr),int(nt)))
    Sinks = np.empty((int(npr),int(nt)))

    for j in range(0,nt,1):
        for i in range(0,npr,1):
            prof = pr[i][j]
            if i == 0:
                p = vg.P1(pset)
            elif i == 1:
                p = vg.P2(pset)
            elif i == 2:
                p = vg.P4(pset)
            else:
                p = vg.P6(pset)
            
            poi = np.flip(prof)
            theta = vg.thetaFun(poi,p)
            MatricPot[:,i,j] = poi
            WaterCon[:,i,j] = theta
            Velocity[:,i,j] = -np.flip(vels[i][j])
            Sinks[i,j] = sink_pr[i][j]
    theta = WaterCon
    ts = 1
    delt = delt/ts
    BigSink = []
    pp = 0.5
    for pro in range(0,4,1):
        if (pro == 0):
            p = vg.P1(pset)
            sy = vg.EstSy(p)
        elif (pro == 1):
            p = vg.P2(pset)
            sy = vg.EstSy(p)
        elif (pro == 2):
            p = vg.P4(pset)
            sy = vg.EstSy(p)
        else:
            p = vg.P6(pset)
            sy = vg.EstSy(p)

        rr = 0
        C_NH4 = np.empty((np.shape(MatricPot)[0],ts*num_time))
        C_NH4[:,0] = cnhic
        C_NO2 = np.empty((np.shape(MatricPot)[0],ts*num_time))
        C_NO2[:,0] = cno2
        C_NO3 = np.empty((np.shape(MatricPot)[0],ts*num_time))
        C_NO3[:,0] = cnic

        lilsink = []
        lilsink.append(a6)
        for j in range(1,nt,1):
            for timestep in range(0,ts,1):
                rr +=1
                for i in range(1,np.shape(MatricPot)[0]-1,1):
                            
                    v = Velocity[i,pro,j-1]/theta[i,pro,j-1]
                    vm1 = Velocity[i-1,pro,j-1]/theta[i-1,pro,j-1]
                    vp1 = Velocity[i+1,pro,j-1]/theta[i+1,pro,j-1]
    
                    # AMMONIUM - NH4
                    # ------------------------------
                    if (i == 1):
                        if ((C_NH4[i+1,rr-1] - C_NH4[i,rr-1]) != 0):
                            roe = (C_NH4[i,rr-1] - C_NH4[i-1,rr-1])/(C_NH4[i+1,rr-1] - C_NH4[i,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    else:
                        if ((C_NH4[i,rr-1] - C_NH4[i-1,rr-1]) != 0):
                            roe = (C_NH4[i-1,rr-1] - C_NH4[i-2,rr-1])/(C_NH4[i,rr-1] - C_NH4[i-1,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    v_avg = 0.5*(v+vm1)
                    Fm = v_avg*theta[i-1,pro,j-1]*C_NH4[i-1,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NH4[i,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i-1,pro,j-1]/2/spac[i]))*theta[i-1,pro,j-1]*C_NH4[i-1,rr-1])
                    if ((C_NH4[i+1,rr-1] - C_NH4[i,rr-1]) != 0):
                        roe = (C_NH4[i,rr-1] - C_NH4[i-1,rr-1])/(C_NH4[i+1,rr-1] - C_NH4[i,rr-1])
                        limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                    else:
                        limiter = 1
                    
                    v_avg = 0.5*(v+vp1)
                    Fp = v_avg*theta[i,pro,j-1]*C_NH4[i,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i+1,pro,j-1]/2/spac[i]))*theta[i+1,pro,j-1]*C_NH4[i+1,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NH4[i,rr-1])
                        
                        
                    no = gaussian(WaterCon[i,pro,j-1]/sy,mu,sig)
                    reactionsNH4 = -no*mu_NH4*C_NH4[i,rr-1]*(C_NH4[i,rr-1]/(K_NH4 + C_NH4[i,rr-1]))
                    retardation = 1 + (rho_b*KOC_NH4/theta[i,pro,j])
                    
                    if (i == 11):
                        ss = pp*Sinks[pro,j-1]*C_NH4[i,rr-1]*sy
                    else:
                        ss = 0
                        
                    if (i<= 10):
                        miner = minerr*gaussian(WaterCon[i,pro,j-1]/sy,mu3,sig3)
                    else:
                        miner = 0
                    
                    C_NH4[i,rr] = (theta[i,pro,j-1]*C_NH4[i,rr-1] + (delt/spac[i])*(Fm-Fp) + delt*(reactionsNH4+miner) - delt*ss)/(retardation*theta[i,pro,j])
                    
                    # NITRITE - NO2
                    # ------------------------------------
                    if (i == 1):
                        if ((C_NO2[i+1,rr-1] - C_NO2[i,rr-1]) != 0):
                            roe = (C_NO2[i,rr-1] - C_NO2[i-1,rr-1])/(C_NO2[i+1,rr-1] - C_NO2[i,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    else:
                        if ((C_NO2[i,rr-1] - C_NO2[i-1,rr-1]) != 0):
                            roe = (C_NO2[i-1,rr-1] - C_NO2[i-2,rr-1])/(C_NO2[i,rr-1] - C_NO2[i-1,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    v_avg = 0.5*(v+vm1)
                    Fm = v_avg*theta[i-1,pro,j-1]*C_NO2[i-1,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NO2[i,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i-1,pro,j-1]/2/spac[i]))*theta[i-1,pro,j-1]*C_NO2[i-1,rr-1])
                    if ((C_NO2[i+1,rr-1] - C_NO2[i,rr-1]) != 0):
                        roe = (C_NO2[i,rr-1] - C_NO2[i-1,rr-1])/(C_NO2[i+1,rr-1] - C_NO2[i,rr-1])
                        limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                    else:
                        limiter = 1
                    
                    v_avg = 0.5*(v+vp1)
                    Fp = v_avg*theta[i,pro,j-1]*C_NO2[i,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i+1,pro,j-1]/2/spac[i]))*theta[i+1,pro,j-1]*C_NO2[i+1,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NO2[i,rr-1])
                
                    no = gaussian(WaterCon[i,pro,j-1]/sy,mu,sig)
                    reactionsNO2 = abs(reactionsNH4) - no*mu_NO2*C_NO2[i,rr-1]*(C_NO2[i,rr-1]/(K_NO2 + C_NO2[i,rr-1]))
                            
                    if (i == 11):
                        ss = pp*Sinks[pro,j-1]*C_NO2[i,rr-1]*sy
                    else:
                        ss = 0
                    
                    C_NO2[i,rr] = (theta[i,pro,j-1]*C_NO2[i,rr-1] + (delt/spac[i])*(Fm-Fp) + delt*reactionsNO2 - delt*ss)/(theta[i,pro,j]) 
            
                    # NITRITE - NO3
                    # ------------------------------------
                    if (i == 1):
                        if ((C_NO3[i+1,rr-1] - C_NO3[i,rr-1]) != 0):
                            roe = (C_NO3[i,rr-1] - C_NO3[i-1,rr-1])/(C_NO3[i+1,rr-1] - C_NO3[i,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    else:
                        if ((C_NO3[i,rr-1] - C_NO3[i-1,rr-1]) != 0):
                            roe = (C_NO3[i-1,rr-1] - C_NO3[i-2,rr-1])/(C_NO3[i,rr-1] - C_NO3[i-1,rr-1])
                            limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                        else:
                            limiter = 1
                    v_avg = 0.5*(v+vm1)
                    Fm = v_avg*theta[i-1,pro,j-1]*C_NO3[i-1,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NO3[i,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i-1,pro,j-1]/2/spac[i]))*theta[i-1,pro,j-1]*C_NO3[i-1,rr-1])
                    if ((C_NO3[i+1,rr-1] - C_NO3[i,rr-1]) != 0):
                        roe = (C_NO3[i,rr-1] - C_NO3[i-1,rr-1])/(C_NO3[i+1,rr-1] - C_NO3[i,rr-1])
                        limiter = max(0.0, min(2*roe,1.0), min(roe,2.0))
                    else:
                        limiter = 1
                    
                    v_avg = 0.5*(v+vp1)
                    Fp = v_avg*theta[i,pro,j-1]*C_NO3[i,rr-1] + \
                        limiter*(((v_avg/2) - (v_avg*delt*Velocity[i+1,pro,j-1]/2/spac[i]))*theta[i+1,pro,j-1]*C_NO3[i+1,rr-1] - \
                                 ((v_avg/2) - (v_avg*delt*Velocity[i,pro,j-1]/2/spac[i]))*theta[i,pro,j-1]*C_NO3[i,rr-1])
                        
                    dent = gaussian(WaterCon[i,pro,j-1]/sy,mu2,sig2)
                    reactionsNO3 = (no*mu_NO2*C_NO2[i,rr-1]*(C_NO2[i,rr-1]/(K_NO2 + C_NO2[i,rr-1]))) - \
                                    dent*mu_NO3*C_NO3[i,rr-1]*(C_NO3[i,rr-1]/(K_NO3 + C_NO3[i,rr-1]))
                            
                    if (i == 11):
                        ss = Sinks[pro,j-1]*C_NO3[i,rr-1]*sy
                        lilsink.append(C_NO3[i,rr-1])
                    else:
                        ss = 0
                        
                    C_NO3[i,rr] = (theta[i,pro,j-1]*C_NO3[i,rr-1] + (delt/spac[i])*(Fm-Fp) + delt*reactionsNO3 - delt*ss)/(theta[i,pro,j]) 
                    
                C_NH4[0,rr] = 2*C_NH4[1,rr] - C_NH4[2,rr]
                C_NO2[0,rr] = 2*C_NO2[1,rr] - C_NO2[2,rr]
                C_NO3[0,rr] = 2*C_NO3[1,rr] - C_NO3[2,rr]
    
                C_NH4[-1,rr] = 2*C_NH4[-2,rr] - C_NH4[-3,rr]
                C_NO2[-1,rr] = 2*C_NO2[-2,rr] - C_NO2[-3,rr]
                C_NO3[-1,rr] = 2*C_NO3[-2,rr] - C_NO3[-3,rr]
                
        BigSink.append(lilsink)
    return BigSink, C_NH4, C_NO3
 
def gaussian(x, muu, sigg):
    return np.exp(-np.power(x - muu, 2.) / (2 * np.power(sigg, 2.)))








