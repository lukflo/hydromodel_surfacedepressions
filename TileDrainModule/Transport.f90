    subroutine Transport()
    
        use global
        implicit none
        
        double precision :: Area, Areap, Aream, v_avg, Fm, Fp, Areanew, ang, roe, limiter
        integer :: drain, ii

        decay = 0.1/86400 ! 1/s 
        E = 0.0527 ! m^2/s - this is a rough upper bound of E for a fully developed turbulent (?) pipe flow. For OC rethink maybe...
        
        do drain = segs,1,-1
            ! FOR NITROGEN            
            if (drain.eq.1) then
                phi(drain,1,2) = 3.185595/Q(drain,1,1)
            else
                phi(drain,1,2) = 0.1312125/Q(drain,1,1)
            endif
            
            ! FOR PHOSPHOROUS
            !phi(drain,1,2) = Conc_BC(drain,days) + 10.0
            
            do ii = 2,(N(drain)-1),1
                if (ii.eq.2) then
                    roe = (phi(drain,ii,1) - phi(drain,ii-1,1))/(phi(drain,ii+1,1) - phi(drain,ii,1))
                    limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                    limiter = max(0.0, min(1.0,roe)) ! minmod
                    if (roe.ne.roe) then
                        limiter =1 
                    end if
                else
                    roe = (phi(drain,ii-1,1) - phi(drain,ii-2,1))/(phi(drain,ii,1) - phi(drain,ii-1,1))
                    limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                    limiter = max(0.0, min(1.0,roe)) ! minmod
                    if (roe.ne.roe) then
                        limiter =1 
                    end if                    
                end if
                if (y(drain,ii,1).gt.D(drain)) then
                    ang = 2*pi
                else
                    ang = 2*acos(1 - 2*y(drain,ii,1)/D(drain))
                end if
                Area = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))  
                
                if (y(drain,ii-1,1).gt.D(drain)) then
                    ang = 2*pi
                else
                    ang = 2*acos(1 - 2*y(drain,ii-1,1)/D(drain))
                end if
                Aream = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))   
                v_avg = 0.5*((Q(drain,ii-1,1)/Aream) + (Q(drain,ii,1)/Area))
                Fm = v_avg*Aream*phi(drain,ii-1,1) + &
                    limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii,1)/Area)/delx(drain)))*Area*phi(drain,ii,1) - &
                    (0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii-1,1)/Aream)/delx(drain)))*Aream*phi(drain,ii-1,1))
                   
                roe = (phi(drain,ii,1) - phi(drain,ii-1,1))/(phi(drain,ii+1,1) - phi(drain,ii,1))
                limiter = max(0.0,min(2*roe,1.0),min(roe,2.0)) ! superbee
                limiter = max(0.0, min(1.0,roe)) ! minmod
                if (roe.ne.roe) then
                    limiter =1 
                end if
                
                if (y(drain,ii+1,1).gt.D(drain)) then
                    ang = 2*pi
                else
                    ang = 2*acos(1 - 2*y(drain,ii+1,1)/D(drain))
                end if
                Areap = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))  
                v_avg = 0.5*((Q(drain,ii+1,1)/Areap) + (Q(drain,ii,1)/Area))
                Fp = v_avg*Area*phi(drain,ii,1) + &
                    limiter*((0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii+1,1)/Areap)/delx(drain)))*Areap*phi(drain,ii+1,1) - &
                    (0.5*v_avg - (0.5*v_avg*delt*(Q(drain,ii,1)/Area)/delx(drain)))*Area*phi(drain,ii,1))              
                
                if (y(drain,ii,1+1).gt.D(drain)) then
                    ang = 2*pi
                else
                    ang = 2*acos(1 - 2*y(drain,ii,1+1)/D(drain))
                end if
                Areanew = (D(drain)*D(drain)/4)*((ang/2) - sin(ang/2)*cos(ang/2))           
                
                phi(drain,ii,2) = (Area*phi(drain,ii,1) + (delt/delx(drain))*(Fm-Fp))/Areanew + &
                    delt*GW_load(drain,ii)/(Areanew*delx(drain)) + delt*ctrib(drain,ii,1)/(Areanew*delx(drain)) 

            end do
            
            phi(drain,N(drain),2) = (Area*phi(drain,N(drain),1) - (delt/delx(drain))*((Q(drain,N(drain),1))*phi(drain,N(drain),1) - &
                (Q(drain,N(drain)-1,1))*phi(drain,N(drain)-1,1)))/Areanew
           

            if (drain.gt.1) then
                ctrib(RecSeg(drain),int(Loca(drain)/delx(drain)),1) = phi(drain,N(drain),2)*Q(drain,N(drain),1) ! units of mg/s     
            end if       
            
            
        end do
    
        deallocate(GW_Load)
        
        
        inquire(file = "last_timestep_C.txt", exist = exist)
        if (exist) then
            open (unit=10200, file="last_timestep_C.txt", status = "replace", action="write")
        else
            open (unit=10200, file="last_timestep_C.txt", status = "new", action="write")
        end if
        do jj = 1,segs
        write(10200,6001) (phi(jj,i,2), i = 1,N(jj))
6001   format ( 5250(16x, e16.8) ) 
        end do
            
            
            
        CLOSE (10200, STATUS='KEEP') 
            
        inquire(file = "OutletC.txt", exist = exist)
        if (exist) then
            open (unit=2010, file="OutletC.txt", status = "old", position = "append", action="write")
        else
            open (unit=2010, file="OutletC.txt", status = "new", action="write")
        end if
    
        write(2010,6003) (phi(1,N(1),2))
6003     format ( 5250(16x, e16.8) )       
            
        inquire(file = "OutlettimeC.txt", exist = exist)
        if (exist) then
            open (unit=2011, file="OutlettimeC.txt", status = "old", position = "append", action="write")
        else
            open (unit=2011, file="OutlettimeC.txt", status = "new", action="write")
        end if
    
        write(2011,6004) (k)
    6004   format ( 5250(16x, e16.8) )  
        
        phi(:,:,1) = phi(:,:,2)
        
    
    end    
        
