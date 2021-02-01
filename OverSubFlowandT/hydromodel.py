import numpy as np
import vanGenuchten as vg
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from landlab.io import read_esri_ascii
from landlab.components.overland_flow import OverlandFlow

def Voyager(j,k,Paths):
    tilenumber = -1
    segnumber = -1
    for p in range(0,len(Paths),1):
        if tilenumber != -1:
            break
        for i in range(0,len(Paths[p]),1):
            if [j,k] == Paths[p][i]:
                tilenumber = p
                segnumber = i
                break
    return [tilenumber, segnumber]

def RichardsModel(psi,t,dz,n,p,vg,qTfun,qBot,psiTop,psiBot,DrainLoc,draincoeff,drain_exp,timestep):
    global savew
    global saveq
    
    C = vg.CFun(psi,p)
    q = np.zeros(n+1)
    w = np.zeros(n+1)

    if psiTop == []:
        if t>(timestep/86400):
            q[n] = qTfun(timestep/86400)
        else:
            q[n] = qTfun(t) 
    else:
        KTop = vg.KFun(np.zeros(1) + psi[n-1], p)
        q[n] = -KTop*((psiTop - psi[n-1])/dz*2 + 1)
#        print(q[n])

    if (psi[DrainLoc] > 0.0):
       w[DrainLoc] = ((draincoeff/vg.EstSy(p))*((psi[DrainLoc])**drain_exp))
#                    
    # Lower boundary
    if qBot == []:
        if psiBot == []:
            # Free drainage
            KBot = vg.KFun(np.zeros(1) + psi[0],p)
            q[0] = -KBot
        else:
            # Type 1 boundary
            KBot = vg.KFun(np.zeros(1) + psiBot,p)
            q[0] = -KBot*((psi[0] - psiBot)/dz*2 + 1.0)    
    else:
        # Type 2 boundary
        q[0] = qBot
        
    # Internal nodes
    i = np.arange(0,n-1)
    Knodes = vg.KFun(psi,p)
    Kmid = (Knodes[i+1] + Knodes[i])/2.0
    
    # Continuity
#    if psiT <= 0:
    j = np.arange(1,n)
    q[j] = -Kmid*((psi[i+1] - psi[i])/dz + 1.0)
    i = np.arange(0,n)
    dpsidt = (-(q[i+1] - q[i])/dz)/C - w[i]/C

#    savew.append(sum(w))
    savew.append(w[DrainLoc])#/C[DrainLoc])
#    savet.append(t)
    saveq.append(q)

    return dpsidt

global savew
global saveq

pset = np.loadtxt('input_param.txt')
ICdepth = pset[-5]
draincoeff = pset[-4]
drain_exp = pset[-3]

DrainLoc = int(round(ICdepth/.1))
dz = 0.1
ProfileDepth = 3
zz = np.arange(dz/2.0, ProfileDepth, dz)
n = zz.size

timestep = 900
t = np.linspace(0,(timestep/86400),2)
qTop = np.zeros(2)
qTfun = interp1d(t,qTop)
qBot = 0 # no flow bottom
psiTop = []
psiBot = []

p = np.loadtxt('rainfall_mm.txt')
precip_md = []
for j in range(0,len(p),1):
    for k in range(0,4,1):
        precip_md.append(p[j]*0.024)


sd1 = 93600 # s
sd2 = 518400
sd3 = 525600
sd4 = 590400
sd5 = 608400
model_run_time = 517500
#model_run_time = 117500
ET = 0.00128/4

PROF = []
GWA = []
TILE = []
VELS = []
for j in range(0,4,1):
    P = np.zeros((n,1))
    po = -zz + ICdepth
    P[:,0]= po
    elapsed_time = 0.
    P1_td = []
    gwa = []
    Prof_1 = []
    vels = []
    tc = 0
    while elapsed_time < model_run_time:
        print(j,elapsed_time)
        save_w = []
        save_q = []
        t = np.linspace(0,(timestep/86400),2) # 10 min time step for one hour
        if j == 0:
            p = vg.P1(pset)
        elif j == 1:
            p = vg.P2(pset)
        elif j == 2:
            p = vg.P4(pset)
        else:
            p = vg.P6(pset)
        
        qTop = np.zeros(2)
        qBot = 0 # no flow bottom
        psiTop = []
        psiBot = []
        GWA_n = -1
        infi = -999
        precip_d = precip_md[tc]
        tc += 1
        if (P[-1,0] < 0): #unsaturated soil column
            if (elapsed_time < sd1):
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                h2 = ProfileDepth + (precip_d*(timestep/86400)) 
                infi = KTop*(h2-(zz[-1]+P[-1,0]))/(dz/2)
                if infi  < precip_d:
                    qTop[0:2] = -infi
                    GWA_n = 1
                else:
                    qTop[0:2] = -precip_d
                    GWA_n = 2
                qTfun = interp1d(t,qTop)
            elif ((elapsed_time > sd2) and (elapsed_time < sd3)):
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                h2 = ProfileDepth + (precip_d*(timestep/86400)) 
                infi = KTop*(h2-(zz[-1]+P[-1,0]))/(dz/2)
    
                if infi  < precip_d:
                    qTop[0:2] = -infi
                    GWA_n = 1
                else:
                    qTop[0:2] = -precip_d
                    GWA_n = 2
                qTfun = interp1d(t,qTop)
                
            elif ((elapsed_time > sd4) and (elapsed_time < sd5)):
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                h2 = ProfileDepth + (precip_d*(timestep/86400)) 
                infi = KTop*(h2-(zz[-1]+P[-1,0]))/(dz/2)
    
                if infi  < precip_d:
                    qTop[0:2] = -infi
                    GWA_n = 1
                else:
                    qTop[0:2] = -precip_d
                    GWA_n = 2
                qTfun = interp1d(t,qTop)
    
            else:
                qTop[0:2]= ET # Evaporation
                qTfun = interp1d(t,qTop) 
                GWA_n = 4
        else: # saturated soil column- runoff
            if (elapsed_time < sd1):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 6
            elif ((elapsed_time > sd2) and (elapsed_time < sd3)):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 6
            elif ((elapsed_time > sd4) and (elapsed_time < sd5)):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    psiTop = []
                    GWA_n = 6
            else:
                qTop[0:2]= ET #  Evaporation
                qTfun = interp1d(t,qTop) 
                GWA_n = 7
                
        savew = []
        saveq = []
        psi0 = P[:,0]
        psi = odeint(RichardsModel,psi0,t,args=(dz,n,p,vg,qTfun,qBot,psiTop,psiBot,DrainLoc,draincoeff,drain_exp,timestep),
                     rtol = 1.49012e-8, atol = 1.49012e-8,mxstep=500000)
    
        if psi[-1][-1] > (dz/2):
            print('Overshot BC')
            if (elapsed_time < sd1):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 6
                   
            elif ((elapsed_time > sd2) and (elapsed_time < sd3)):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 6
            elif ((elapsed_time > sd4) and (elapsed_time < sd5)):
                psiTop = 0
                KTop=vg.KFun(np.zeros(1)+P[-1,0],p)[0]
                PressFlux = KTop*((psiTop - P[-1,0])/dz*2 + 1)
                if PressFlux > precip_d:
                    qTop[0:2] = -precip_d
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 5
                else:
                    qTop[0:2] = -PressFlux
                    qTfun = interp1d(t,qTop) 
                    GWA_n = 6
            else:
                qTop[0:2]= ET # Evaporation
                qTfun = interp1d(t,qTop) 
                GWA_n = 7
    
            savew = []
            saveq = []   
            psi0 = P[:,0]                             
            psi = odeint(RichardsModel,psi0,t,args=(dz,n,p,vg,qTfun,qBot,psiTop,psiBot,DrainLoc,draincoeff,drain_exp,timestep),
                         rtol = 1.49012e-10, atol = 1.49012e-10, mxstep=500000)
        P[:,0] = psi[-1] 
        
        P1_td.append(savew[-1]/86400)
        vels.append(saveq[-1])
        gwa.append([GWA_n,infi])
        Prof_1.append(psi[-1])
        elapsed_time+= timestep
        
    PROF.append(Prof_1)
    GWA.append(gwa)
    TILE.append(P1_td)
    VELS.append(vels)

raise Exception

##########################################################################
# Soils
soils = np.load('SOILArray.npy')
intakes = np.load('intakes2.npy')
intakes_ = np.load('intakes_nu.npy')
intakes2 = []
for j in range(0,len(intakes_),1):
    if j != 1:
        intakes2.append(intakes_[j])

ocoeff = pset[-2]
    
wa = []
intake_index = []

for i in range(0,len(intakes),2):
    ppp = 0
    for j in range(np.shape(soils)[0]-1,-1,-1):
        for k in range(0,np.shape(soils)[1],1):
            if ((j == (72-intakes[i][0])) and (k == intakes[i][1])):
                intake_index.append(ppp+62)
                break
            ppp +=1
    
p = np.loadtxt('rainfall_mm.txt')
precip_md = []
for j in range(0,len(p),1):
    for k in range(0,4,1):
        precip_md.append(p[j]*0.024)

doi1 = pset[-1]
doi1 = .07

T_h = []
cc = 0
# OV Flow grid
(mg,z) = read_esri_ascii('gstc240_2.txt', name ='topographic__elevation')
mg.set_watershed_boundary_condition(z, -9999)
of = OverlandFlow(mg, mannings_n = 0.025, alpha = 0.4, steep_slopes = True,h_init = 1e-7)

elapsed_time = 0.
ET2 = -0.0006/86400

H = np.zeros((4464,205000)) 
qq = np.zeros((8794,205000)) 
hm = 0

H[:,0] = mg.at_node['surface_water__depth']
qq[:,0] = mg.at_link['surface_water__discharge']

tt = []
tc = 0
while elapsed_time < model_run_time:
    water_intake = []
    ppp = 0
    precip_d = precip_md[tc]
    tc += 1
    t = np.linspace(0,(timestep/86400),2) # 10 min time step for one hour
    if (elapsed_time < sd1):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    elif ((elapsed_time > sd2) and (elapsed_time < sd3)):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    elif ((elapsed_time > sd4) and (elapsed_time < sd5)):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    else:
        gridded_water_array = (ET2)*np.ones(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    for j in range(np.shape(soils)[0]-1,-1,-1):
        for k in range(0,np.shape(soils)[1],1):
            if z[ppp] > 0:
                if ((soils[j,k]==1)or(soils[j,k]==2)or(soils[j,k]==4)or(soils[j,k]==7)or(soils[j,k]==13)or(soils[j,k]==14)or(soils[j,k]==17)):
                    p = vg.P1(pset)
                    tt.append(1)
                elif ((soils[j,k]==8)or(soils[j,k]==11)or(soils[j,k]==12)or(soils[j,k]==16)or(soils[j,k]==20)):
                    p = vg.P2(pset)
                    tt.append(2)
                elif (soils[j,k]==15):
                    p = vg.P4(pset)
                    tt.append(4)
                else:
                    p = vg.P6(pset)
                    tt.append(6)

                if tt[-1] == 1:
                    GWA_n = GWA[0][int(elapsed_time/timestep)][0]
                    infi = GWA[0][int(elapsed_time/timestep)][1]
                elif tt[-1] == 2:
                    GWA_n = GWA[1][int(elapsed_time/timestep)][0]
                    infi = GWA[1][int(elapsed_time/timestep)][1]
                elif tt[-1] == 4:
                    GWA_n = GWA[2][int(elapsed_time/timestep)][0]
                    infi = GWA[2][int(elapsed_time/timestep)][1]
                else:
                    GWA_n = GWA[3][int(elapsed_time/timestep)][0]
                    infi = GWA[3][int(elapsed_time/timestep)][1]
                    
                print(tt[-1],elapsed_time, ppp,GWA_n)

                if (GWA_n == 1):
                    if elapsed_time > 1800:
                        excess = (precip_d - infi)/(86400)
                        gridded_water_array[ppp] += excess
                        hm+= excess
                    else:
                        gridded_water_array[ppp] += 0
                elif (GWA_n == 2): # 
                    gridded_water_array[ppp] += 0
                elif (GWA_n == 3):
                    psiTop = 0
                    KTop=vg.KFun(np.zeros(1)+P[-1,ppp],p)[0]                                
                    gridded_water_array[ppp] = -KTop*((psiTop - psi0[n-1])/dz*2 + 1)/86400
                elif (GWA_n == 4):
                    gridded_water_array[ppp] += 0
                elif (GWA_n == 5):
                    gridded_water_array[ppp] += 0
                elif (GWA_n == 6):
                    gridded_water_array[ppp] += starting_precip
                    hm += starting_precip
                elif (GWA_n == 7):
                    gridded_water_array[ppp] += 0
                else:
                    raise Exception
                
                if (ppp in intake_index):
                    if (mg.at_node['surface_water__depth'][ppp] > doi1):
                        gridded_water_array[ppp] += -((ocoeff*((2*9.81*(mg.at_node['surface_water__depth'][ppp]-doi1))**4.05))/(90*90))
                        water_intake.append([ppp,mg.at_node['surface_water__depth'][ppp],((ocoeff*((2*9.81*(mg.at_node['surface_water__depth'][ppp]-doi1))**4.05)))])
                    else:
                        gridded_water_array[ppp] += 0
                        water_intake.append([ppp, mg.at_node['surface_water__depth'][ppp], 0])
            ppp = ppp + 1 
    ov_elap = 0
    while (ov_elap < timestep):        
        cc = cc + 1
        mg.at_node['surface_water__depth'] += gridded_water_array*5 
        
        for j in range(0,len(mg.at_node['surface_water__depth']),1):
            if mg.at_node['surface_water__depth'][j] < 0:
                mg.at_node['surface_water__depth'][j] = 1e-7
                
        of.run_one_step(dt = 5)
        ov_elap += 5
        H[:,cc] = mg.at_node['surface_water__depth']
        qq[:,cc] = mg.at_link['surface_water__discharge'] 
        T_h.append(ov_elap + elapsed_time)
  
        print(ov_elap)

    elapsed_time +=  timestep
    wa.append(water_intake)
