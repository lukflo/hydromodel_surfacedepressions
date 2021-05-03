!  TileDrainSRC.f90 
!
!  FUNCTIONS:
!  TileDrainSRC - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: TileDrainSRC
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program TileDrainSRC

    use global 
    !USE OMP_LIB
    implicit none

    open(100, file = 'Input\Param.txt')
    open(101, file = 'Input\BC_rev.txt')
    open(1022, file = 'Input\lq.txt')
    open(1023, file = 'Input\ly.txt')
    open(601, file = 'Input\wt.txt')
    
    read(100,*)
    read(100,*) segs, stressperiods, delt, drain_coeff, exp_coeff, depth_to_drain, transport_flag
    read(100,*)
    
    delt1 = delt
    
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


    allocate(Q(segs,maxval(N),2))
    allocate(y(segs,maxval(N),2))
    allocate(qtrib(segs,maxval(N),1))
    allocate(Qlast(segs,maxval(N)))
    allocate(ylast(segs,maxval(N)))
    Q(:,:,:) = -99
    y(:,:,:) = -99
    qtrib(:,:,:) = 0 

   ! ------------------------------------------------
2894 continue 
     
    allocate(Q_up(segs,stressperiods))
    allocate(y_down(segs,stressperiods))
    allocate(IC_Q(segs,maxval(N)))
    allocate(IC_y(segs,maxval(N)))
    do jj = 1,segs,1
        read(1022,*) IC_Q(jj,1:N(jj))
        read(1023,*) IC_y(jj,1:N(jj))
        read(101,*) Q_up(jj,:)
    end do    
    
    Q(:,:,1) = IC_Q(:,:)
    y(:,:,1) = IC_y(:,:)
    Qlast(:,:) = Q(:,:,1)
    ylast(:,:) = y(:,:,1)

    do jj = segs,2,-1
        qtrib(int(RecSeg(jj)),int(Loca(jj)/delx(jj)),1) = -Q(jj,N(jj),1)/delx(int(RecSeg(jj)))
    end do

    last_trib = -1
    ! ------------------------------------------------
    do days = 1,stressperiods
        
        allocate(wt_day(segs,maxval(N)))
        do jj = 1,segs,1
            read(601,*) wt_day(jj,:)
        end do
    
        k= 0.0
        loop_count = -1
        
        do while (k.lt.iii)
            
            loop_count = loop_count + 1

            print*,'SP: ', days, 'TS: ', k, 'LC: ', loop_count
            
            tol_trib = tol_trib1 
            tol_y = tol_y1
            tol_Q = tol_Q1
            
            do jj= 1,segs
                tol_trib(jj) = abs(0.01*Q(jj,N(jj),1)) ! use to be 0.005
            end do
            
            ! Courant-Freidrich-Lewy condition            
            delt1 = (6000.0/segs)*delx(1)/(Q(1,maxval(N)-5,1)/((D(1)*D(1)/4)*(acos(1 - 2*y(1,maxval(N)-5,1)/D(1)) -&
                sin(acos(1 - 2*y(1,maxval(N)-5,1)/D(1)))*cos(acos(1 - 2*y(1,maxval(N)-5,1)/D(1))))))
            delt = delt1
            tol_y(:) = 0.01
            tol_Q(:) = 0.01

2893        continue             
            
            converge = segs
            iter_count = 0
            do while (converge.gt.0)
                
                if (iter_count.gt.30) then
                    if (delt.gt.(0.01)) then
                        delt = 0.3*delt
                        solv_count = 0
                        iter_count = 0
                        do jj = segs,2,-1
                            qtrib(int(RecSeg(jj)),int(Loca(jj)/delx(jj)),1) = -Q(jj,N(jj),1)/delx(int(RecSeg(jj)))
                        end do                          
                        print*, 'ITER COUNTER TIME ADJUST ', last_trib, delt
                        go to 2893
                    else
                        tol_trib(last_trib) = 5*tol_trib(last_trib)
                        iter_count = 0
                        solv_count = 0
                        do jj = segs,2,-1
                             qtrib(int(RecSeg(jj)),int(Loca(jj)/delx(jj)),1) = -Q(jj,N(jj),1)/delx(int(RecSeg(jj)))
                        end do                          
                        print*, 'ADJUSTING TRIB_TOL FOR TRIB ', last_trib
                        go to 2893
                    end if
                end if

                allocate(wt(N(converge)))
                allocate(reach_off(2*N(converge))) 
                allocate(guessQ(N(converge))) 
                allocate(guessy(N(converge))) 
                allocate(oldQ(N(converge))) 
                allocate(oldy(N(converge)))
                allocate(J(2*N(converge),2*N(converge)))
                allocate(Je(2*N(converge),4))
                allocate(resid(2*N(converge))) 
                allocate(x(2*N(converge)))

                ! read in drain flow ql - dim(1)
                wt(:) = wt_day(converge,:)
              
                if (k.eq.(0.0)) then
                    guessQ = Q(converge,:,1) - Qlast(converge,:) + Q(converge,:,1)
                    guessy = y(converge,:,1) - ylast(converge,:) + y(converge,:,1)
                else
                    guessQ = ((Q(converge,:,1) - Qlast(converge,:))/olddelt)*delt + Q(converge,:,1)
                    guessy = ((y(converge,:,1) - ylast(converge,:))/olddelt)*delt + y(converge,:,1)
                end if
         
                
                ! previous time step
                oldQ = Q(converge,:,1)
                oldy = y(converge,:,1)
                ! Adjust DS BC for tributaries
                if (converge.gt.1) then
                    y_down(converge,days) = y(int(RecSeg(converge)),int(Loca(converge)/delx(converge)),2)
                end if  
            
                ! Create Jacobian        
                call MatrixAssembly()
                call JacobianInversion()
                
                flag = 1
                solv_count = 0
                do while (flag == 1)
                    ! check which reaches are within error tolerance
                    do r = 1,(2*N(converge)),1
                        
                        if (mod(r,2).ne.0) then
                            if (abs(x(r)).gt.tol_Q(converge)*abs(oldQ(int((r/2)+1)))) then
                                reach_off(r) = 1
                            else
                                reach_off(r) = 0
                            end if
                        else          
                            if (abs(x(r)).gt.tol_y(converge)*abs(oldy(int((r/2))))) then
                                reach_off(r) = 1
                            else
                                reach_off(r) = 0
                            end if
                        end if
                    end do 
                    
                    if (sum(reach_off).eq.0) then
                        do r = 1,(2*N(converge)),1
                            if (mod(r,2).ne.0) then
                                guessQ((r/2)+1) = guessQ((r/2)+1) + x(r)
                            else
                                guessy((r/2)) = guessy((r/2)) + x(r)
                                if (guessy((r/2)).lt.(0.00000000001)) then
                                    guessy((r/2)) = oldy((r/2))
                                end if
                            end if
                        end do
                        flag = 0
                        solv_count = 0
                    else
                        do r = 1,(2*N(converge)),1
                            if (mod(r,2).ne.0) then
                                guessQ((r/2)+1) = guessQ((r/2)+1) + x(r)
                            else
                                guessy((r/2)) = guessy((r/2)) + x(r)
                                if (guessy((r/2)).lt.(0.00000000001)) then
                                    guessy((r/2)) = oldy((r/2))
                                end if
                            end if
                        end do
                        
                        call MatrixAssembly()
                        call JacobianInversion()
                                            
                        solv_count = solv_count + 1
                end if
                    
                end do                     
                if (converge.lt.(segs+1)) then
                    Q(converge,:,2) = guessQ
                    y(converge,:,2) = guessy
                    converge = converge - 1
                else
                    if (abs(guessQ(N(converge)) - -1*qtrib(int(RecSeg(converge)),int(Loca(converge)/delx(converge)),1)* &
                        delx(int(RecSeg(converge)))).lt.tol_trib(converge)) then
                        Q(converge,:,2) = guessQ
                        y(converge,:,2) = guessy
                        converge = converge - 1
                    else          
                        qtrib(int(RecSeg(converge)),int(Loca(converge)/delx(converge)),1) = &
                            -guessQ(N(converge))/delx(int(RecSeg(converge)))
                        last_trib = converge
                        converge = segs
                        iter_count = iter_count + 1
                    end if 
                end if
            
                deallocate(reach_off)
                deallocate(guessQ)
                deallocate(guessy)
                deallocate(oldQ)
                deallocate(oldy)
                deallocate(J)
                deallocate(Je)
                deallocate(resid)
                deallocate(x)
                deallocate(wt)
                    
            end do                
            
            Qlast(:,:) = Q(:,:,1)
            ylast(:,:) = y(:,:,1)
            
            y(:,:,1) = y(:,:,2)
            Q(:,:,1) = Q(:,:,2)
            k = k + delt
            
            do jj = segs,2,-1
                qtrib(int(RecSeg(jj)),int(Loca(jj)/delx(jj)),1) = -Q(jj,N(jj),1)/delx(int(RecSeg(jj)))
            end do
            
            olddelt = delt   
                       
            inquire(file = "OutletQ2.txt", exist = exist)
            if (exist) then
                open (unit=201, file="OutletQ2.txt", status = "old", position = "append", action="write")
            else
                open (unit=201, file="OutletQ2.txt", status = "new", action="write")
            end if
    
            write(201,60003) (Q(1,N(1),2))
60003       format ( 5250(16x, e16.8) )       
            
            inquire(file = "Outlettime2.txt", exist = exist)
            if (exist) then
                open (unit=201, file="Outlettime2.txt", status = "old", position = "append", action="write")
            else
                open (unit=201, file="Outlettime2.txt", status = "new", action="write")
            end if
    
            write(201,60004) (k)
    60004   format ( 5250(16x, e16.8) )      
            
            
        end do
        
       deallocate(wt_day)
        
        
                inquire(file = "last_day_Qnosi.txt", exist = exist)
        if (exist) then
            open (unit=1022, file="last_day_Qnosi.txt", status = "old", position = "append", action="write")
        else
            open (unit=1022, file="last_day_Qnosi.txt", status = "new", action="write")
        end if
        do jj = 1,segs
        write(1022,61000) (Q(jj,i,2), i = 1,N(jj))
61000   format ( 5250(16x, e16.8) ) 
        end do
            
            
        inquire(file = "last_day_ynosi.txt", exist = exist)
        if (exist) then
            open (unit=1023, file="last_day_ynosi.txt", status = "old", position = "append", action="write")
        else
            open (unit=1023, file="last_day_ynosi.txt", status = "new", action="write")
        end if
        do jj = 1,segs
        write(1023,61001) (y(jj,i,2), i = 1,N(jj))
61001   format ( 5250(16x, e16.8) ) 
        end do          
            
        CLOSE (1022, STATUS='KEEP') 
        CLOSE (1023, STATUS='KEEP') 
        CLOSE (231, STATUS='KEEP') 
            
    end do
  
    end program TileDrainSRC

  
