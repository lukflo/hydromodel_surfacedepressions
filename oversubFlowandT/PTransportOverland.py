import numpy as np

def SedTranCap(qs,n,h,kt,gamma,bet,dx):   
    V = qs/dx
    Sf = abs((V*n/(h**(2/3)))**2)
    Je = abs(kt*(qs**bet))*(Sf**gamma)
    return Je

manning = 0.4
bet = 2.2
gam = 0.8
kt = 5857
soils = np.load('SOILArray.npy')
# hardcoded, below
diff_coeff = 0
nx =  62
xmax = 90*nx
ny = 72
ymax = 90*ny
dx = xmax/(nx)
dy = ymax/(ny)
#
#q = np.load('_qq2nu.npy')
H = np.load('_Hhsinu.npy')
si = np.load('_Enssinu.npy')
nt = 50000

intakes = np.load('intakes2.npy')
in_p = []
for k in range(0,400,1):
    in_p.append(si[k])
intake_index = []
sinkarray = np.zeros((ny,nx,nt))
for i in range(0,len(intakes),2):
    ppp = 0
    for j in range(np.shape(soils)[0]-1,-1,-1):
        for k in range(0,np.shape(soils)[1],1):
            if ((j == (72-intakes[i][0])) and (k == intakes[i][1])):
                intake_index.append([72-j,k,ppp+62])
                break
            ppp +=1
mult = -1
for j in range(0,nt,1):
    if np.mod(j,180) == 0:
        mult += 1
    for k in range(0,len(intake_index),1):
        sinkarray[intake_index[k][0],intake_index[k][1],j] = in_p[mult][k][2]

R = np.load('Right2.npy')
L = np.load('Left2.npy')
U = np.load('Up2.npy')
D = np.load('Down2.npy')
raise Exception
c = np.zeros((ny,nx,nt)) 
bigSink = []
nn = -1
elapsed_time = 0
for n in range(0,30000,2):
    print(n,nt)
    elapsed_time+=10
    delt = 10
    ii = 0
    nn += 1
    for j in range(0,ny,1):
        for i in range(0,nx,1):
            if ((i == 0) or (j == 0) or (i == (nx-1)) or (j == (ny-1))):
                ii += 1
            else:
                v_x = 0.5*(R[nn][j,i]+L[nn][j,i])/90
                v_x_p1 = 0.5*(R[nn][j,i+1]+L[nn][j,i+1])/90
                v_x_m1 = 0.5*(R[nn][j,i-1]+L[nn][j,i-1])/90
                v_y = 0.5*(U[nn][j,i]+D[nn][j,i])/90
                v_y_p1 = 0.5*(U[nn][j+1,i]+D[nn][j+1,i])/90
                v_y_m1 = 0.5*(U[nn][j-1,i]+D[nn][j-1,i])/90
                
                if (i == 1):
                    if ((c[j,i+1,nn] - c[j,i,nn]) != 0):
                        roe = (c[j,i,nn] - c[j,i-1,nn])/(c[j,i+1,nn] - c[j,i,nn])
                        limiterx = min(0.0,min(2*roe,1),min(roe,2.0))
                    else:
                        limiterx = 1
                else:
                    if ((c[j,i,nn] - c[j-1,i,nn]) != 0):
                        roe = (c[j-1,i,nn] - c[j-2,i,nn])/(c[j,i,nn] - c[j-1,i,nn])
                        limiterx = min(0.0,min(2*roe,1),min(roe,2.0))
                    else:
                        limiterx = 1
                        
                if (j == 1):
                    if ((c[j,i+1,nn] - c[j,i,nn]) != 0):
                        roe = (c[j,i,nn] - c[j,i-1,nn])/(c[j,i+1,nn] - c[j,i,nn])
                        limitery = min(0.0,min(2*roe,1),min(roe,2.0))
                    else:
                        limitery = 1
                else:
                    if ((c[j,i,nn] - c[j,i-1,nn]) != 0):
                        roe = (c[j,i-1,nn] - c[j,i-2,nn])/(c[j,i,nn] - c[j,i-1,nn])
                        limitery = min(0.0,min(2*roe,1),min(roe,2.0))
                    else:
                        limitery = 1
                        
                v_avgx = 0.5*(v_x + v_x_m1)
                v_avgy = 0.5*(v_y + v_y_m1)
                
                Fmx = v_avgx*c[j-1,i,nn] + limiterx*(((v_avgx/2) - (v_avgx*delt*v_x/2/dx))*c[j,i,nn] - \
                              ((v_avgx/2) - (v_avgx*delt*v_x_m1/2/dx))*c[j-1,i,nn])
                Fmy = v_avgy*c[j,i-1,nn] + limitery*(((v_avgy/2) - (v_avgy*delt*v_y/2/dy))*c[j,i,nn] - \
                              ((v_avgy/2) - (v_avgy*delt*v_y_m1/2/dy))*c[j,i-1,nn]) 
                
                if ((c[j+1,i,nn] - c[j,i,nn]) != 0):
                    roe = (c[j,i,nn] - c[j-1,i,nn])/(c[j+1,i,nn] - c[j,i,nn])
                    limiterx = min(0.0,min(2*roe,1),min(roe,2.0))
                else:
                    limiterx = 1
    
                if ((c[j,i+1,nn] - c[j,i,nn]) != 0):
                    roe = (c[j,i,nn] - c[j,i-1,nn])/(c[j,i+1,nn] - c[j,i,nn])
                    limitery = min(0.0,min(2*roe,1),min(roe,2.0))
                else:
                    limitery = 1
                    
                v_avgx = 0.5*(v_x + v_x_p1)
                v_avgy = 0.5*(v_y + v_y_p1)
                
                diffusiony = (-(2*diff_coeff*delt/(dy*dy)))*c[j,i,nn] +\
                    (diff_coeff*delt/(dy*dy))*(c[j+1,i,nn] + c[j-1,i,nn])

                diffusionx = (-(2*diff_coeff*delt/(dx*dx)))*c[j,i,nn] + \
                    (diff_coeff*delt/(dx*dx))*(c[j,i+1,nn] + c[j,i-1,nn])

                Fpx = v_avgx*c[j,i,nn] + limiterx*(((v_avgx/2) - (v_avgx*delt*v_x_p1/2/dx))*c[j+1,i,nn] - \
                              ((v_avgx/2) - (v_avgx*delt*v_x/2/dx))*c[j,i,nn])
                Fpy = v_avgy*c[j,i,nn] + limitery*(((v_avgy/2) - (v_avgy*delt*v_y_p1/2/dy))*c[j,i+1,nn] - \
                              ((v_avgy/2) - (v_avgy*delt*v_y/2/dy))*c[j,i,nn]) 
                
                ripq = np.sqrt((v_x*dx)**2 + (v_y*dy)**2)

                if (elapsed_time < 95400):
                    c[j,i,nn+1] = c[j,i,nn] + (delt/dx)*(Fmx - Fpx) + (delt/dy)*(Fmy - Fpy) + \
                                delt*SedTranCap(ripq,manning,H[ii,n],kt,gam,bet,dx) - \
                                delt*sinkarray[j,i,nn]*c[j,i,nn+1]/dx/dy + diffusionx + diffusiony
                elif ((elapsed_time > 518400) and (elapsed_time < 525600)):
                    kt = 457142
                    c[j,i,nn+1] = c[j,i,nn] + (delt/dx)*(Fmx - Fpx) + (delt/dy)*(Fmy - Fpy) + \
                                delt*SedTranCap(ripq,manning,H[ii,n],kt,gam,bet,dx) - \
                                delt*sinkarray[j,i,nn]*c[j,i,nn+1]/dx/dy + diffusionx + diffusiony
                elif ((elapsed_time > 590400) and (elapsed_time < 608400)):
                    kt = 457152
                    c[j,i,nn+1] = c[j,i,nn] + (delt/dx)*(Fmx - Fpx) + (delt/dy)*(Fmy - Fpy) + \
                                delt*SedTranCap(ripq,manning,H[ii,n],kt,gam,bet,dx) - \
                                delt*sinkarray[j,i,nn]*c[j,i,nn+1]/dx/dy + diffusionx + diffusiony
                else: # bad way to account for particule settling...
                    kt = 457142
                    c[j,i,nn+1] = c[j,i,nn] + (delt/dx)*(Fmx - Fpx) + (delt/dy)*(Fmy - Fpy) + \
                                    delt*SedTranCap(ripq,manning,H[ii,n],kt,gam,bet,dx) - \
                                    delt*sinkarray[j,i,nn]*c[j,i,nn+1]/dx/dy - \
                                    c[j,i,nn]*.001 + diffusionx + diffusiony

                if sinkarray[j,i,n] > 0:
                    bigSink.append([j,i,int(n),c[j,i,nn+1]])
                ii += 1

inc_v = np.zeros((42,nt))
for t in range(0,np.shape(bigSink)[0],1):
    print(t)
    for k in range(0,np.shape(intake_index)[0],1):
        if ((bigSink[t][0] == intake_index[k][0]) and (bigSink[t][1] == intake_index[k][1])):
            inc_v[k,int(bigSink[t][2])] = bigSink[t][3] #mg/l
            break   

I_c = np.zeros((42,16200))
for j in range(0,np.shape(inc_v)[0],1):
    print(j)
    for k in range(0,np.shape(inc_v)[1],1):
        if (np.mod(k,180) == 0):
            I_c[j,int(k/180)] = abs(inc_v[j,k])
       
            
                
