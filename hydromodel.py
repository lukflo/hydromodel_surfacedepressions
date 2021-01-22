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
       w[DrainLoc] = (((1/7)*draincoeff/vg.EstSy(p))*((psi[DrainLoc])**drain_exp))
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
ICdepth = 1.0386*pset[-5]
draincoeff = 1.239*pset[-4]
drain_exp = 2.437*pset[-3]

ICdepth = 1.0546*pset[-5] # meant for 2
draincoeff = 1.397*pset[-4]
drain_exp = 1.916*pset[-3]

#ICdepth = 1.017*pset[-5] # meant for 4
#draincoeff = 1.2547*pset[-4]
#drain_exp = 1.737*pset[-3]

ICdepth = 1.00601*pset[-5] # meant for 4
draincoeff = 1.112*pset[-4]
drain_exp = 1.51658*pset[-3]
#
#ICdepth = 1.0544*pset[-5]  #meant for 2
#draincoeff = 1.4135*pset[-4]
#drain_exp = 1.89*pset[-3]

DrainLoc = int(round(ICdepth/.1))
dz = 0.1
ProfileDepth = 3
zz = np.arange(dz/2.0, ProfileDepth, dz)
n = zz.size
#
#np1 = 1777
#np2 = 544
#np4 = 60
#np6 = 10

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

PROF = []
GWA = []
TILE = []
VELS = []
for j in range(0,4,1):
    P = np.zeros((n,1))
    po = -zz + ICdepth
    P[:,0]= po
    elapsed_time = 0.
    storm_duration = 93600 # s
    model_run_time = 972000# s
    model_run_time = 517500
    precip_d = 0.067
    timestep = 900
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
            if (elapsed_time < storm_duration):
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
            elif ((elapsed_time > 518400) and (elapsed_time < 525600)):
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
                
            elif ((elapsed_time > 590400) and (elapsed_time < 608400)):
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
                qTop[0:2]= 0.00128/4 # Evaporation
                qTfun = interp1d(t,qTop) 
                GWA_n = 4
        else: # saturated soil column- runoff
            if (elapsed_time < storm_duration):
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
            elif ((elapsed_time > 518400) and (elapsed_time < 525600)):
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
            elif ((elapsed_time > 590400) and (elapsed_time < 608400)):
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
                qTop[0:2]= 0.00128/4 #  Evaporation
                qTfun = interp1d(t,qTop) 
                GWA_n = 7
                
        savew = []
        saveq = []
        psi0 = P[:,0]
        psi = odeint(RichardsModel,psi0,t,args=(dz,n,p,vg,qTfun,qBot,psiTop,psiBot,DrainLoc,draincoeff,drain_exp,timestep),
                     rtol = 1.49012e-8, atol = 1.49012e-8,mxstep=500000)
    
        if psi[-1][-1] > 0.05:
            print('Overshot BC')
            if (elapsed_time < storm_duration):
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
                   
            elif ((elapsed_time > 518400) and (elapsed_time < 525600)):
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
            elif ((elapsed_time > 590400) and (elapsed_time < 608400)):
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
                qTop[0:2]= 0.00127/4 # Evaporation
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


##########################################################################
# try with 4
# 1.09, 1.585, 1.012    
    
#np1 = 1110
#np2 = 648
#np4 = 328
#np6 = 305
#np1 = 1777
#np2 = 544
#np4 = 61
#np6 = 10
#S = []
#for k in range(0,450,2):
#    S1 = 90*90*TILE[0][k]
#    
#    S2 = 90*90*TILE[1][k]
#    
#    S4 = 90*90*TILE[2][k]
#    
#    S6 = 90*90*TILE[3][k]
#
#    S.append((S1*np1)+(S2*np2)+(S4*np4)+(S6*np6))
#     
#S = [i/30 for i in S]
#import matplotlib.pyplot as plt
#plt.plot(S)
#plt.grid()
#raise Exception

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
storm_duration = 93600 # s
model_run_time = 972000 # s
model_run_time = 517500
#starting_precip = (2.7778e-7)*2.79 # m/s
#precip_d = 0.067
timestep = 900

# Initialize (probably a better way to do this...)
#H = np.zeros((4464,205000)) 
#qq = np.zeros((8794,205000)) 
hm = 0

#H[:,0] = mg.at_node['surface_water__depth']
#qq[:,0] = mg.at_link['surface_water__discharge']

tt = []
tc = 0
while elapsed_time < model_run_time:
    water_intake = []
    ppp = 0
    precip_d = precip_md[tc]
    tc += 1
    t = np.linspace(0,(timestep/86400),2) # 10 min time step for one hour
    if (elapsed_time < 95400):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    elif ((elapsed_time > 518400) and (elapsed_time < 525600)):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    elif ((elapsed_time > 590400) and (elapsed_time < 608400)):
        gridded_water_array = np.zeros(np.shape(mg.at_node['surface_water__depth']))
        starting_precip = precip_d/86400
    else:
        gridded_water_array = (-0.0006/86400)*np.ones(np.shape(mg.at_node['surface_water__depth']))
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
                        gridded_water_array[ppp] += -((0.1*0.12*0.0015*ocoeff*((2*9.81*(mg.at_node['surface_water__depth'][ppp]-doi1))**3.9))/(90*90))
                        water_intake.append([ppp,mg.at_node['surface_water__depth'][ppp],((0.1*0.12*0.0015*ocoeff*((2*9.81*(mg.at_node['surface_water__depth'][ppp]-doi1))**3.9)))])
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
#                    print('sdsdf')
                mg.at_node['surface_water__depth'][j] = 1e-7
                
        of.run_one_step(dt = 5)
#        of.overland_flow(dt = 5)
        ov_elap += 5
    
#        H[:,cc] = mg.at_node['surface_water__depth']
        T_h.append(ov_elap + elapsed_time)
#        qq[:,cc] = mg.at_link['surface_water__discharge'] 
  
        print(ov_elap)

    elapsed_time +=  timestep
    wa.append(water_intake)
#np.save('_Hhsinu.npy',H)
#np.save('_qq2nu.npy',qq)
#np.save('_Enssinu.npy',water_intake)   
import matplotlib.pyplot as plt

#si = []
#for j in range(0,375,1):
#    print(j)
#    s = 0
#    for k in range(0,42,1):
#        s += wa[j][k][2]
#    si.append(s)
#
##
#plt.plot(si)
#raise Exception
##import matplotlib.pyplot as plt
##from matplotlib.colors import LogNorm
##base = np.load('_Hhsinu.npy')
##BigGrid_h1 = []
###BigGrid_h2 = []
##for f in range(0,90001,1):
##    print(f)
##    smallgrid_h1 = np.zeros((72,62))
###    smallgrid_h2 = np.zeros((72,62))
##    node = 0
##    for j in range(0,np.shape(smallgrid_h1)[0],1):
##        for k in range(0,np.shape(smallgrid_h1)[1],1):
##            smallgrid_h1[j,k] = 100*base[node,f]
###            smallgrid_h2[j,k] = 100*ex1[node,f]
##            node = node + 1
##    BigGrid_h1.append(smallgrid_h1) 
###    BigGrid_h2.append(smallgrid_h2) 
##
##for a in range(0,110000,90):
##    BigGrid1 = np.ma.masked_where(BigGrid_h1[a] ==  BigGrid_h1[a][0,0], BigGrid_h1[a])
##    im=plt.imshow(BigGrid1, norm= LogNorm(vmin = 1e-5, vmax=5))
##    plt.title(a*5)
##    plt.gca().invert_yaxis()
##    plt.pause(0.05)
##    del(im)


Paths = np.load('Paths_Long.npy',allow_pickle = True)
Qtimes = np.zeros((480*len(Paths),len(Paths[0])))
BCs = np.zeros((len(Paths),480))
ows = np.load('outofWS.npy')
for times in range(0,480,1): # 1080
    Ql = -99*np.ones((len(Paths),len(Paths[0])))
    for j in range(1,len(Paths),1):
        Ql[j,len(Paths[j])::] = 0
    print(times)
    inta = wa[times]
    ppp = 0
    tt = []
    ppp_r = 0 
    Psy = []
    hou = []
    bd_tile_time = []
    for j in range(np.shape(soils)[0]-1,-1,-1):
        for k in range(0,np.shape(soils)[1],1):
            if z[ppp] > 0:
                hou.append(Voyager(j,k,Paths))
                for pdumb in range(0,len(inta),1):
                    if inta[pdumb][0] == ppp:
                        if hou[-1][0] != -1:
                            bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                            break
                        else:
                            if (Voyager(j+1,k,Paths)[0] != -1):
                                hou.append(Voyager(j+1,k,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                                break
                            elif (Voyager(j,k+1,Paths)[0] != -1):
                                hou.append(Voyager(j,k+1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                                break
                            elif (Voyager(j-1,k,Paths)[0] != -1):
                                hou.append(Voyager(j-1,k,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                                break
                            elif (Voyager(j,k-1,Paths)[0] != -1):
                                hou.append(Voyager(j,k-1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]]) 
                                break
                            elif (Voyager(j+1,k-1,Paths)[0] != -1):
                                hou.append(Voyager(j+1,k-1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                                break
                            elif (Voyager(j-1,k-1,Paths)[0] != -1):
                                hou.append(Voyager(j-1,k-1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])
                                break
                            elif (Voyager(j+1,k+1,Paths)[0] != -1):
                                hou.append(Voyager(j+1,k+1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]])  
                                break
                            elif (Voyager(j-1,k+1,Paths)[0] != -1):
                                hou.append(Voyager(j-1,k+1,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]]) 
                                break
                            elif (Voyager(j+2,k,Paths)[0] != -1):
                                hou.append(Voyager(j+2,k,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]]) 
                                break 
                            elif (Voyager(j,k+2,Paths)[0] != -1):
                                hou.append(Voyager(j,k+2,Paths))
                                bd_tile_time.append([hou[-1][0],inta[pdumb][2]]) 
                                break 
                            else:
                                print(j,k)
                                raise Exception
##
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
#                    
#                if ((soils[j,k]==1)or(soils[j,k]==5)or(soils[j,k]==11)or(soils[j,k]==12)or(soils[j,k]>=15)):
#                    p = vg.P1(pset)
#                    tt.append(1)
#                elif ((soils[j,k]==2)or(soils[j,k]==6)or(soils[j,k]==10)or(soils[j,k]==13)or(soils[j,k]==14)):
#                    p  = vg.P2(pset)
#                    tt.append(2)
#                elif (soils[j,k]==3):
#                    p = vg.P3()
#                    tt.append(3)
#                elif (soils[j,k]==7):
#                    p = vg.P4(pset)
#                    tt.append(4)
#                elif (soils[j,k]==8):
#                    p = vg.P5()
#                    tt.append(5)
#                elif (soils[j,k]==9):
#                    p = vg.P6(pset)
#                    tt.append(6)
#                else:
#                    p = vg.P2(pset)
#                    tt.append(2)    
#                    
#                    
                if tt[-1] == 1:
                    p = TILE[0][times]
                elif tt[-1] == 2:
                    p = TILE[1][times]
                elif tt[-1] == 4:
                    p = TILE[2][times]
                else:
                    p = TILE[3][times]

                if (hou[-1] != [-1,-1]):
                    Ql[hou[-1][0],hou[-1][1]] = 3*3*p

                    
                Psy.append(p)
                ppp_r += 1
            ppp += 1
        

    for j in range(0,len(ows),1):
        Ql[ows[j][0],ows[j][1]] = Ql[ows[j][0]-1,ows[j][1]]

    for jj in range(0,np.shape(Ql)[0],1):
        for kk in range(0,np.shape(Ql)[1],1):
            if Ql[jj,kk] == -99:
                Ql[jj,kk] = 0
        
    Qtimes[times*len(Paths):((times+1)*len(Paths)),:] = Ql
    # PROBLEM SPECIFIC
    for b in range(0,len(bd_tile_time),1):
        if ((bd_tile_time[b][0] == 276)):
            BCs[239,times] += 1.0*bd_tile_time[b][1]
#            BCs[240,times] += 0.25*bd_tile_time[b][1]
        elif (bd_tile_time[b][0] == 277):
            BCs[248,times] += 0.5*bd_tile_time[b][1]
            BCs[263,times] += 0.5*bd_tile_time[b][1]
        else:
            BCs[bd_tile_time[b][0],times] += abs(bd_tile_time[b][1])

QT = np.zeros((480*len(Paths),30*len(Paths[0])))
for j in range(0,np.shape(Qtimes)[0],1):
    print(j)
    for k in range(0,np.shape(Qtimes)[1],1):
        QT[j,(k*30):((k*30)+30)] = Qtimes[j][k]

for k in range(0,480,1):
    if BCs[0,k] < 0.0002123730:
        BCs[0,k] = 0.0002123730
        
for j in range(0,335,1):        
    for k in range(0,480,1):
        if BCs[j,k] < 0.0000087475:
            BCs[j,k] = 0.0000087475
    
np.savetxt('wt_.txt',QT,fmt='%1.10f')
np.savetxt('BC_.txt',BCs,fmt='%1.10f')



    
