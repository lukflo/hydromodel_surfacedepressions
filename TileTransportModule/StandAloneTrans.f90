!  StandAloneTrans.f90 
!
!  FUNCTIONS:
!  StandAloneTrans - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: StandAloneTrans
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program StandAloneTrans

    implicit none

    ! Variables
    double precision, dimension(:), allocatable :: delx, theta, L, D, mann, So, tol_y, Loca, tol_trib, tol_trib1
    double precision, dimension(:), allocatable :: tol_y1, tol_Q, tol_Q1
    integer, dimension(:), allocatable :: chantype, RecSeg, DSBCType, USBCType, N
    
    double precision :: Area, Areap, Aream, v_avg, Fm, Fp, Areanew, ang, roe, limiter, delt, delt_t
    double precision :: drain_coeff, exp_coeff, depth_to_drain, transport_flag, k, best_t
    double precision :: slopey, slopeq, y_smooth_nui, y_smooth_i, y_smooth_pi, y_smooth_mi
    double precision :: Q_smooth_nui, Q_smooth_i, Q_smooth_pi, Q_smooth_mi
    integer :: drain, ii, stressperiods, days,jj, segs, loop_count, i, iii
    real :: g = 9.81, pi = 3.14159265359, mult_num
    
    double precision, dimension(:,:,:), allocatable :: phi, phi3
    double precision, dimension(:,:), allocatable :: wt_day, GW_Conc, GW_Load, Conc_BC, Q, y, ctrib, Qn, yn
    double precision, dimension(:,:), allocatable :: Q3, y3, Qn3, yn3, wt_day3, GW_Conc3, GW_load3

    logical :: exist

    open(100,file='Input\Param.txt')
    open(601, file = 'Input\wt_.txt')
    open(602, file = 'Input\conc.txt')
    open(701, file = 'last_day_Q.txt')
    open(702, file = 'last_day_y.txt')
    open(1022, file = 'Input\ICc2.txt')

    read(100,*)
    read(100,*) segs, stressperiods, delt, drain_coeff, exp_coeff, depth_to_drain, transport_flag
    read(100,*)

    allocate(delx(segs))
    allocate(theta(segs))
    allocate(L(segs))
    allocate(chantype(segs))
    allocate(D(segs))
    allocate(mann(segs))
    allocate(So(segs))
    allocate(tol_y(segs))
    allocate(tol_y1(segs))
    allocate(tol_Q(segs))
    allocate(tol_Q1(segs))
    allocate(RecSeg(segs))
    allocate(DSBCType(segs))
    allocate(Loca(segs))
    allocate(USBCType(segs))
    allocate(N(Segs))
    allocate(tol_trib(Segs))
    allocate(tol_trib1(Segs))
    
    do jj = 1,segs,1
        read(100,*) delx(jj), theta(jj), L(jj), chantype(jj), D(jj), mann(jj), &
            So(jj), tol_y(jj), RecSeg(jj), DSBCType(jj), Loca(jj), USBCType(jj),&
            tol_trib(jj), tol_Q(jj)
    end do  
    tol_trib1 = tol_trib
    tol_y1 = tol_y
    tol_Q1 = tol_Q
    do jj = 1,segs,1
        N(jj) = L(jj)/delx(jj)
    end do

    delx = delx/3
    
    allocate(phi(segs,maxval(N),2))
    allocate(phi3(segs,3*maxval(N),2))
    allocate(ctrib(segs,3*maxval(N)))
    allocate(Conc_BC(segs,stressperiods)) 
    
    ! THIS SECTION IS NEEDED FOR NITROGEN
    !################################
    do jj = 1,segs,1
        read(1022,*) phi3(jj,1:(3*N(jj)),1)
    end do   
    !###################################
    
    ctrib(:,:) = 0

    allocate(Q(segs,maxval(N)))
    allocate(y(segs,maxval(N)))
    allocate(Q3(segs,3*maxval(N)))
    allocate(y3(segs,3*maxval(N)))
    
    do jj = 1,segs,1
        read(701,*) Q(jj,1:N(jj))
        read(702,*) y(jj,1:N(jj))
    end do

    do jj = 1,segs,1
        do ii = 1,maxval(N),1
            Q3(jj,3*(ii-1)+1) = Q(jj,ii) 
            Q3(jj,3*(ii-1)+2) = Q(jj,ii) 
            Q3(jj,3*(ii-1)+3) = Q(jj,ii) 
            y3(jj,3*(ii-1)+1) = y(jj,ii) 
            y3(jj,3*(ii-1)+2) = y(jj,ii) 
            y3(jj,3*(ii-1)+3) = y(jj,ii) 
        end do
    end do
    
    do days = 1,stressperiods

        allocate(wt_day(segs,maxval(N)))
        allocate(GW_Conc(segs,maxval(N))) ! UNITS OF MG/m^3
        !allocate(GW_Load(segs,maxval(N)))
        allocate(wt_day3(segs,3*maxval(N)))
        allocate(GW_Conc3(segs,3*maxval(N))) ! UNITS OF MG/m^3
        allocate(GW_Load3(segs,3*maxval(N)))
        !
        do jj = 1,segs,1
            read(601,*) wt_day(jj,:)
            read(602,*) GW_conc(jj,:)
        end do
        
        do jj = 1,segs,1
            do ii = 1,maxval(N),1
                wt_day3(jj,3*(ii-1)+1) = wt_day(jj,ii) 
                wt_day3(jj,3*(ii-1)+2) = wt_day(jj,ii) 
                wt_day3(jj,3*(ii-1)+3) = wt_day(jj,ii) 
                GW_conc3(jj,3*(ii-1)+1) = GW_conc(jj,ii) 
                GW_conc3(jj,3*(ii-1)+2) = GW_conc(jj,ii) 
                GW_conc3(jj,3*(ii-1)+3) = GW_conc(jj,ii) 
            end do
        end do
        
        allocate(Qn(segs,maxval(N)))
        allocate(yn(segs,maxval(N)))
        allocate(Qn3(segs,3*maxval(N)))
        allocate(yn3(segs,3*maxval(N)))

        do jj = 1,segs,1
            read(701,*) Qn(jj,1:N(jj))
            read(702,*) yn(jj,1:N(jj))
        end do
        
        do jj = 1,segs,1
            do ii = 1,maxval(N),1
                Qn3(jj,3*(ii-1)+1) = Qn(jj,ii) 
                Qn3(jj,3*(ii-1)+2) = Qn(jj,ii) 
                Qn3(jj,3*(ii-1)+3) = Qn(jj,ii) 
                yn3(jj,3*(ii-1)+1) = yn(jj,ii) 
                yn3(jj,3*(ii-1)+2) = yn(jj,ii) 
                yn3(jj,3*(ii-1)+3) = yn(jj,ii) 
            end do
        end do
        
        do jj = 1,segs,1
            do ii = 1,maxval(N),1
        !        GW_load(jj,ii) = 1000*wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj)
                GW_load3(jj,3*(ii-1)+1) = wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj)
        !        GW_load3(jj,3*(ii-1)+1) = 0
        !        !GW_load3(jj,3*(ii-1)+1) = max(0.0,0.05*1000*wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj))
                GW_load3(jj,3*(ii-1)+2) = wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj)
        !        GW_load3(jj,3*(ii-1)+2) = 0
        !        !GW_load3(jj,3*(ii-1)+1) = max(0.0,0.05*0.1000*wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj))
                GW_load3(jj,3*(ii-1)+3) = wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj)
        !        GW_load3(jj,3*(ii-1)+3) = 0
        !        !GW_load3(jj,3*(ii-1)+1) = max(0.0,0.05*1000*wt_day(jj,ii)*(GW_conc(jj,ii))*delx(jj))
            end do
        end do
        
        mult_num = 0.99
        
2893 continue
        
        delt = 20.0
        do drain = segs,1,-1   
            do ii = 2,((3*N(drain))-1),1
                if (days.lt.80) then
                    delt_t = abs(mult_num*delx(drain)/(abs(Qn3(drain,ii))/(D(drain)*yn3(drain,ii))))
                else
                    delt_t = abs(0.97*mult_num*delx(drain)/(abs(Qn3(drain,ii))/(D(drain)*yn3(drain,ii))))
                end if
                if (delt_t.lt.delt) then
                    delt = delt_t
                end if
            end do
        end do

        k= 0.0
        loop_count = -1
        iii = 900.0
        do while (k.lt.iii)
            k = k + delt
            loop_count = loop_count + 1

            print*,'SP: ', days, 'TS: ', k, 'LC: ', loop_count
            

            do drain = segs,1,-1                    
            !do drain = 1,1,-1
                 !FOR NITROGEN            
                if (drain.eq.1) then
                    phi3(drain,1,2) = ((15.0/15.0)*3.185595/Qn3(drain,1))
                else
                    phi3(drain,1,2) = ((15.0/15.0)*0.1312125/Qn3(drain,1))
                endif

                ! FOR PHOSPHOROUS                
                !phi3(drain,1,2) = 50.0 + max(1000.0*Conc_BC(drain,days),0.0)


                do ii = 2,((3*N(drain))-1),1
                    
                    slopey = (yn3(drain,ii) - y3(drain,ii))/900.0
                    y_smooth_nui = slopey*k + y3(drain,ii)
                    slopey = (yn3(drain,ii) - y3(drain,ii))/900.0
                    y_smooth_i = slopey*(k-delt) + y3(drain,ii)
                    slopey = (yn3(drain,ii+1) - y3(drain,ii+1))/900.0
                    y_smooth_pi = slopey*(k-delt) + y3(drain,ii+1)
                    slopey = (yn3(drain,ii-1) - y3(drain,ii-1))/900.0
                    y_smooth_mi = slopey*(k-delt) + y3(drain,ii-1)

                    slopeq = (Qn3(drain,ii) - Q3(drain,ii))/900.0
                    q_smooth_nui = slopeq*k + Q3(drain,ii)
                    slopeq = (Qn3(drain,ii) - Q3(drain,ii))/900.0
                    q_smooth_i = slopeq*(k-delt) + Q3(drain,ii)
                    slopeq = (Qn3(drain,ii+1) - Q3(drain,ii+1))/900.0
                    q_smooth_pi = slopeq*(k-delt) + Q3(drain,ii+1)
                    slopeq = (Qn3(drain,ii-1) - Q3(drain,ii-1))/900.0
                    q_smooth_mi = slopeq*(k-delt) + Q3(drain,ii-1)     
                    
                    !if (((q_smooth_i/(D(drain)*y_smooth_i))*delt/delx(drain)).gt.0.9) then
                    !    print*, drain, ii, (q_smooth_i/(D(drain)*y_smooth_i))*delt/delx(drain)
                    !    pause
                    !end if
                    
                    if (ii.eq.2) then
                        roe = (phi3(drain,ii,1) - phi3(drain,ii-1,1))/(phi3(drain,ii+1,1) - phi3(drain,ii,1))
                        limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                        limiter = max(0.0, min(1.0,roe)) ! minmod
                        if (roe.ne.roe) then
                            limiter =1 
                        end if
                    else
                        roe = (phi3(drain,ii-1,1) - phi3(drain,ii-2,1))/(phi3(drain,ii,1) - phi3(drain,ii-1,1))
                        limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                        limiter = max(0.0, min(1.0,roe)) ! minmod
                        if (roe.ne.roe) then
                            limiter =1 
                        end if                    
                    end if
                    
                    !limiter = 1
                    
                    !if (y(drain,ii).gt.D(drain)) then
                    !    ang = 2*pi
                    !else
                    !    ang = 2*acos(1 - 2*y(drain,ii)/D(drain))
                    !end if
                    !Area = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))  
                    !Area = D(drain)*y(drain,ii)
                    Area = D(drain)*y_smooth_i
                
                    !if (y(drain,ii-1).gt.D(drain)) then
                    !    ang = 2*pi
                    !else
                    !    ang = 2*acos(1 - 2*y(drain,ii-1)/D(drain))
                    !end if
                    !Aream = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))   
                    !Aream = D(drain)*y(drain,ii-1)
                    Aream = D(drain)*y_smooth_mi
                    
                    !v_avg = 0.5*((Q(drain,ii-1)/Aream) + (Q(drain,ii)/Area))
                    !Fm = v_avg*Aream*phi(drain,ii-1,1) + &
                    !    limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii)/Area)/delx(drain)))*Area*phi(drain,ii,1) - &
                    !    (0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii-1)/Aream)/delx(drain)))*Aream*phi(drain,ii-1,1))
                    
                    v_avg = 0.5*((Q_smooth_mi/Aream) + (Q_smooth_i/Area))
                    Fm = v_avg*Aream*phi3(drain,ii-1,1) + &
                        limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q_smooth_i/Area)/delx(drain)))*Area*phi3(drain,ii,1) - &
                        (0.5*v_avg - (0.5*v_avg*delt*(Q_smooth_mi/Aream)/delx(drain)))*Aream*phi3(drain,ii-1,1))
                   
                    roe = (phi3(drain,ii,1) - phi3(drain,ii-1,1))/(phi3(drain,ii+1,1) - phi3(drain,ii,1))
                    limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                    limiter = max(0.0, min(1.0,roe)) ! minmod
                    if (roe.ne.roe) then
                        limiter =1 
                    end if
                    
                    !limiter = 1
                
                    !if (y(drain,ii+1).gt.D(drain)) then
                    !    ang = 2*pi
                    !else
                    !    ang = 2*acos(1 - 2*y(drain,ii+1)/D(drain))
                    !end if
                    !Areap = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))  
                    !Areap = D(drain)*y(drain,ii+1)
                    Areap = D(drain)*y_smooth_pi
                    
                    !v_avg = 0.5*((Q(drain,ii+1)/Areap) + (Q(drain,ii)/Area))
                    !Fp = v_avg*Area*phi(drain,ii,1) + &
                    !    limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii+1)/Areap)/delx(drain)))*Areap*phi(drain,ii+1,1) - &
                    !    (0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii)/Area)/delx(drain)))*Area*phi(drain,ii,1))   
                    
                    v_avg = 0.5*((Q_smooth_pi/Areap) + (Q_smooth_i/Area))
                    Fp = v_avg*Area*phi3(drain,ii,1) + &
                        limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q_smooth_pi/Areap)/delx(drain)))*Areap*phi3(drain,ii+1,1) - &
                        (0.5*v_avg - (0.5*v_avg*delt*(Q_smooth_i/Area)/delx(drain)))*Area*phi3(drain,ii,1))   
                
                    !if (y(drain,ii).gt.D(drain)) then
                    !    ang = 2*pi
                    !else
                    !    ang = 2*acos(1 - 2*y(drain,ii)/D(drain))
                    !end if
                    !Areanew = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))     
                    !Areanew = D(drain)*y(drain,ii)
                    Areanew = D(drain)*y_smooth_nui

                    phi3(drain,ii,2) = (Area*phi3(drain,ii,1) + (delt/delx(drain))*(Fm-Fp))/Areanew + &
                        delt*GW_load3(drain,ii)/(Areanew*delx(drain)) + delt*ctrib(drain,ii)/(Areanew*delx(drain)) 
                    
                    if (phi3(drain,ii,2).ne.phi3(drain,ii,2)) then
                        mult_num = 0.05*mult_num
                        print*, 'INSTABILITY',drain, delt
                        pause
                        go to 2893
                    end if

                    
                end do    
                
                phi3(drain,3*N(drain),2) = 2*phi3(drain,3*(N(drain))-1,2) - phi3(drain,(3*N(drain))-2,2)
           
                !if (phi3(drain,3*N(drain),2).lt.0.0) then
                !    print*, drain, ii
                !    pause
                !    phi3(drain,3*N(drain),2) = 0
                !end if
                
                
                if (drain.gt.1) then
                    !ctrib(RecSeg(drain),int(Loca(drain)/delx(drain))) = phi(drain,N(drain),2)*Q(drain,N(drain)) ! units of mg/s
                    ctrib(RecSeg(drain),int(Loca(drain)/delx(drain))) = phi3(drain,3*N(drain),2)*Q_smooth_pi! units of mg/s   
                end if 
            end do
            
            
        !    inquire(file = "last_timestep_C.txt", exist = exist)
        !    if (exist) then
        !        open (unit=10200, file="last_timestep_C.txt", status = "replace", action="write")
        !    else
        !        open (unit=10200, file="last_timestep_C.txt", status = "new", action="write")
        !    end if
        !    do jj = 1,segs
        !    write(10200,6001) (phi3(jj,i,2), i = 1,3*N(jj))
        !    !write(10200,6001) (GW_load(1,i), i = 1,N(jj))
        !6001   format ( 5250(16x, e16.8) ) 
        !    end do
        !    
        !    CLOSE (10200, STATUS='KEEP') 
            
            inquire(file = "OutletCP.txt", exist = exist)
            if (exist) then
                open (unit=2010, file="OutletCP.txt", status = "old", position = "append", action="write")
            else
                open (unit=2010, file="OutletCP.txt", status = "new", action="write")
            end if
    
            write(2010,6003) (phi3(1,3*N(1),2))
        6003     format ( 5250(16x, e16.8) )       
            
            inquire(file = "OutlettimeCP.txt", exist = exist)
            if (exist) then
                open (unit=2011, file="OutlettimeCP.txt", status = "old", position = "append", action="write")
            else
                open (unit=2011, file="OutlettimeCP.txt", status = "new", action="write")
            end if
    
            write(2011,6004) (k)
        6004   format ( 5250(16x, e16.8) )  
        
            phi3(:,:,1) = phi3(:,:,2)
            
        end do
        
!        inquire(file = "last_day_c.txt", exist = exist)
!        if (exist) then
!            open (unit=1022, file="last_day_c.txt", status = "old", position = "append", action="write")
!        else
!            open (unit=1022, file="last_day_c.txt", status = "new", action="write")
!        end if
!        write(1022,61000) (phi3(3,i,2), i = 1,N(3))
!61000   format ( 5250(16x, e16.8) ) 
!        
!                inquire(file = "last_day_q.txt", exist = exist)
!        if (exist) then
!            open (unit=1023, file="last_day_q.txt", status = "old", position = "append", action="write")
!        else
!            open (unit=1023, file="last_day_q.txt", status = "new", action="write")
!        end if
!        write(1023,62000) (Q(3,i), i = 1,N(3))
!62000   format ( 5250(16x, e16.8) ) 
!        
!                inquire(file = "last_day_y.txt", exist = exist)
!        if (exist) then
!            open (unit=1024, file="last_day_y.txt", status = "old", position = "append", action="write")
!        else
!            open (unit=1024, file="last_day_y.txt", status = "new", action="write")
!        end if
!        write(1024,63000) (y(3,i), i = 1,N(3))
!63000   format ( 5250(16x, e16.8) ) 
        
!                inquire(file = "last_day_gw.txt", exist = exist)
!        if (exist) then
!            open (unit=1025, file="last_day_gw.txt", status = "old", position = "append", action="write")
!        else
!            open (unit=1025, file="last_day_gw.txt", status = "new", action="write")
!        end if
!        write(1025,64000) (gw_load(3,i), i = 1,N(3))
!64000   format ( 5250(16x, e16.8) ) 
        
        deallocate(wt_day)
        deallocate(GW_Conc)
        !deallocate(GW_Load)
        
        Q = Qn
        y = yn
        
        deallocate(Qn)
        deallocate(yn) 
        
        deallocate(wt_day3)
        deallocate(GW_Conc3)
        deallocate(GW_Load3)
        
        Q3 = Qn3
        y3 = yn3
        
        deallocate(Qn3)
        deallocate(yn3) 
        
    end do
        

            
            
    end program StandAloneTrans
    
    