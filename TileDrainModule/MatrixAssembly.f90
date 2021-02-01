    subroutine MatrixAssembly()
    
        use global
        implicit none
            
        double precision :: Mass_i, Mo_i, B1, Bn
        double precision :: Mass_yi, Mass_yi2, Mass_Qi, Mass_Qi2
        double precision :: Mo_yi, Mo_yi2, Mo_Qi, Mo_Qi2
        double precision :: B1_y1, B1_Q1, Bn_yn, Bn_Qn
        double precision :: A_i_guess, A_i2_guess, A_i_old, A_i2_old
        double precision :: cen_i_guess, cen_i2_guess, cen_i_old, cen_i2_old
        double precision :: Rh_i_guess, Rh_i2_guess, Rh_i_old, Rh_i2_old
        double precision :: A_hi, A_hi2, Rh_hi, Rh_hi2, cen_hi, cen_hi2
        
        double precision :: Sf_yi, Sf_yi2, Sf_Qi, Sf_Qi2
        double precision :: Sf_i_new, Sf_i_old, Sf_i2_new, Sf_i2_old, y_i, So_i
        integer :: ii, col
        double precision :: c, Tw, slopem, Qup_smooth
        double precision :: T1, T2, T3
        
        c = 2*delt/delx(converge)     
        col = 1
        slot_width = 0.01*D(converge)
        
        allocate(ql(N(converge)))
        ql(:) = 0            
        
        !if (Q_upbase(converge,104).eq.Q_upbase(converge,2)) then
        !    slopem = (min(Q_up(converge,days+1), Q_upbase(converge,days+1)) - &
        !        min(Q_up(converge,days),Q_upbase(converge,days)))/900.0
        !    Qup_smooth = slopem*k + min(Q_up(converge,days), Q_upbase(converge,days))
        !    B1 = guessQ(1) - Qup_smooth
        !else
            slopem =( Q_up(converge,days+1) -Q_up(converge,days))/900.0
            Qup_smooth = slopem*k + Q_up(converge,days)
            B1 = guessQ(1) - Qup_smooth
        !end if

        !B1 = guessQ(1)- Q_upbase(converge,1)

        
        if (converge.ge.1) then !for now always assume critical depth DSBC for Seg1
            if (chantype(converge).gt.0) then !OC
                Bn = guessQ(N(converge)) - sqrt(g*((guessy(N(converge))*D(converge))**3)/D(converge))
                Bn_Qn = 1
                Bn_yn = sqrt(g*((guessy(N(converge))*D(converge))**3)/D(converge))*(-1.5/guessy(N(converge)))
            else ! pipe
                if (guessy(N(converge)).gt.D(converge)) then ! pressure
                    Bn = guessQ(N(converge)) - &
                        sqrt(g/slot_width)*(pi*((D(converge)/2)**2) + slot_width*(guessy(N(converge)) - D(converge)))**(3./2.)
                    Bn_Qn = 1
                    Bn_yn = sqrt(g/slot_width)*(pi*((D(converge)/2)**2) + slot_width*(guessy(N(converge)) - D(converge)))**(3./2.)*&
                        ((-1.5*slot_width/(pi*((D(converge)/2)**2) + slot_width*(guessy(N(converge)) - D(converge)))) + &
                        (1./(2*slot_width)))
                else ! open drain
                    alpha = 2*acos(1 - (2*guessy(N(converge))/D(converge)))
                    Bn = guessQ(N(converge)) - &
                        sqrt(g/(D(converge)*sin(alpha/2)))*(((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2)))**(3./2.)
                    Bn_Qn = 1
                    Bn_yn = sqrt(g/(D(converge)*sin(alpha/2)))*(((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2)))**(3./2.)*&
                        ((-1.5*D(converge)*sin(alpha/2)/(((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2)))) + &
                        (1./(2*D(converge)*sin(alpha/2)))*0.5*D(converge)*cos(alpha/2)*(4/D(converge)/sqrt(-(2/D(converge))*guessy(N(converge))*&
                        ((2/D(converge))*guessy(N(converge)) - 2))))
                end if
            end if   
        else ! depth bc
            Bn = guessy(N(converge)) - y_down(converge,days)
            Bn_yn = 1
            Bn_Qn = 0
        end if
        
        ! US Derivatives
        B1_y1 = 0
        B1_Q1 = 1
        
        J(:,:) = 0
        J(1,1) = B1_Q1
        J(1,2) = B1_y1
        J(2*N(converge),(2*N(converge))-1) = Bn_Qn
        J(2*N(converge),2*N(converge)) = Bn_yn
        
        ! Edited Jacobian
        Je(:,:) = 0
        Je(1,3) = B1_Q1
        Je(1,4) = B1_y1
        Je(2*N(converge),1) = Bn_Qn
        Je(2*N(converge),2) = Bn_yn

        resid(1) = -B1
        resid(2*N(converge)) = -Bn

        do ii = 1,(N(converge)-1),1
            
            !if ((converge.eq.291).and.(ii.lt.200)) then
            !    theta(converge) = 0.9
            !end if
        
            if (chantype(converge).gt.0) then ! OC            
                A_i_guess = D(converge)*guessy(ii)
                A_i2_guess = D(converge)*guessy(ii+1)
                A_i_old = D(converge)*oldy(ii)
                A_i2_old = D(converge)*oldy(ii+1)
            
                Rh_i_guess = A_i_guess/(D(converge) + 2*guessy(ii))
                Rh_i2_guess = A_i2_guess/(D(converge) + 2*guessy(ii+1))
                Rh_i_old = A_i_old/(D(converge) + 2*oldy(ii))
                Rh_i2_old = A_i2_old/(D(converge) + 2*oldy(ii+1))
                
                cen_i_guess = 0.5*guessy(ii)*guessy(ii)*D(converge)
                cen_i2_guess = 0.5*guessy(ii+1)*guessy(ii+1)*D(converge)
                cen_i_old = 0.5*oldy(ii)*oldy(ii)*D(converge)
                cen_i2_old = 0.5*oldy(ii+1)*oldy(ii+1)*D(converge)
    
                A_hi = D(converge)
                A_hi2 = D(converge)
                Rh_hi = (D(converge)**2)/((D(converge)+2*guessy(ii))**2)
                Rh_hi2 = (D(converge)**2)/((D(converge)+2*guessy(ii+1))**2)
                cen_hi = guessy(ii)*D(converge)
                cen_hi2 = guessy(ii+1)*D(converge)
                ql(ii) = wt(ii)
                ql(ii+1) = wt(ii+1)
                
            else ! pipe
                if (guessy(ii).ge.D(converge)) then ! pressurized
                    A_i_guess = pi*((D(converge)/2)**2) + slot_width*(guessy(ii) - D(converge))
                    Rh_i_guess = A_i_guess/(pi*D(converge) + 2*(guessy(ii)-D(converge)))
                    A_hi = slot_width
                    Rh_hi = ((pi*D(converge) + 2*(guessy(ii)-D(converge)))*A_hi - 2*A_i_guess)/((pi*D(converge) + &
                        2*(guessy(ii)-D(converge)))**2)
                    
                    !cen_i_guess = 
                    !cen_hi = 
            
                    ql(ii) = wt(ii)
!                    ql(ii) = 0
                    !GW_load(converge,ii) = GW_conc(converge,ii)*wt(ii)*delx(ii)
                else                    
                    alpha = 2*acos(1 - (2*guessy(ii)/D(converge)))
                    A_i_guess = ((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2))
                    Rh_i_guess = A_i_guess/(D(converge)*alpha/2)
                    A_hi = ((D(converge)**2)/8)*(1-cos(alpha))*((4/D(converge))/sqrt((-2*guessy(ii)/D(converge))*&
                        ((2*guessy(ii)/D(converge)) - 2)))
                    
                    !cen_i_guess = 
                    !cen_hi = 
                    
                    Rh_hi = ((D(converge)*alpha/2)*A_hi - A_i_guess*(D(converge)/2)*((4/D(converge))/&
                        sqrt((-2*guessy(ii)/D(converge))*((2*guessy(ii)/D(converge)) - 2))))/((D(converge)*alpha/2)**2)

                    ql(ii) = wt(ii)
!                    ql(ii) = 0
                    !GW_load(converge,ii) = GW_conc(converge,ii)*wt(ii)*delx(ii)
                end if
                
                if (guessy(ii+1).ge.D(converge)) then
                    A_i2_guess = pi*((D(converge)/2)**2) + slot_width*(guessy(ii+1) - D(converge))
                    Rh_i2_guess = A_i2_guess/(pi*D(converge) + 2*(guessy(ii+1)-D(converge)))
                    A_hi2 = slot_width
                    Rh_hi2 = ((pi*D(converge) + 2*(guessy(ii+1)-D(converge)))*A_hi2 - 2*A_i2_guess)/&
                        ((pi*D(converge) + 2*(guessy(ii+1)-D(converge)))**2)
                    
                    !cen_i2_guess = 
                    !cen_hi2 = 
                    
                    ql(ii+1) = wt(ii+1)
                    !ql(ii+1) = 0
                    !GW_load(converge,ii+1) = GW_conc(converge,ii+1)*wt(ii+1)*delx(ii+1)
                else
                    alpha = 2*acos(1 - (2*guessy(ii+1)/D(converge)))
                    A_i2_guess = ((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2))
                    Rh_i2_guess = A_i2_guess/(D(converge)*alpha/2)
                    A_hi2 = ((D(converge)**2)/8)*(1-cos(alpha))*((4/D(converge))/&
                        sqrt((-2*guessy(ii+1)/D(converge))*((2*guessy(ii+1)/D(converge)) - 2)))
                    Rh_hi2 = ((D(converge)*alpha/2)*A_hi2 - A_i2_guess*(D(converge)/2)*((4/D(converge))/&
                        sqrt((-2*guessy(ii+1)/D(converge))*((2*guessy(ii+1)/D(converge)) - &
                        2))))/((D(converge)*alpha/2)**2)
                    
                    !cen_i2_guess = 
                    !cen_hi2 = 
                    
                    ql(ii+1) = wt(ii+1)
                    !ql(ii+1) = 0
                   ! GW_load(converge,ii+1) = GW_conc(converge,ii+1)*wt(ii+1)*delx(ii+1)
                end if
                    
                if (oldy(ii).gt.D(converge)) then
                    A_i_old = pi*((D(converge)/2)**2) + slot_width*(oldy(ii) - D(converge))
                    Rh_i_old = A_i_old/(pi*D(converge) + 2*(oldy(ii)-D(converge)))
                    !cen_i_old = 
                else
                    alpha = 2*acos(1 - (2*oldy(ii)/D(converge)))
                    A_i_old = ((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2))
                    Rh_i_old = A_i_old/(D(converge)*alpha/2)
                    !cen_i_old = 
                end if
                  
                if (oldy(ii+1).gt.D(converge)) then
                    !cen_i2_old = 
                    A_i2_old = pi*((D(converge)/2)**2) + slot_width*(oldy(ii+1) - D(converge))
                    Rh_i2_old = A_i2_old/(pi*D(converge) + 2*(oldy(ii+1)-D(converge)))
                else
                    alpha = 2*acos(1 - (2*oldy(ii+1)/D(converge)))
                    A_i2_old = ((D(converge)**2)/4)*((alpha/2) - sin(alpha/2)*cos(alpha/2))
                    Rh_i2_old = A_i2_old/(D(converge)*alpha/2)
                    !cen_i2_old = 
                end if
            end if
            
                
            ! Co. of Mass Function Evalutation
            Mass_i = A_i_guess + A_i2_guess - A_i_old - A_i2_old + &
                c*(theta(converge)*(guessQ(ii+1)-guessQ(ii)) + (1-theta(converge))*(oldQ(ii+1) - oldQ(ii))) + &
                delt*(theta(converge)*(qtrib(converge,ii,1) + qtrib(converge,ii+1,1)) +&
                (1-theta(converge))*(qtrib(converge,ii,1) + qtrib(converge,ii+1,1))) - &
                delt*(theta(converge)*(ql(ii) + ql(ii+1)) + (1-theta(converge))*(ql(ii) + ql(ii+1)))
            
            ! Co. of Mass Derivatives
            Mass_yi = A_hi
            Mass_yi2 = A_hi2
            Mass_Qi = -c*theta(converge)
            Mass_Qi2 = c*theta(converge)

            
            ! Co. of  Momentum Function Evaluation
            Sf_i_new = ((mann(converge)**2)*guessQ(ii)*abs(guessQ(ii))/(A_i_guess*A_i_guess*Rh_i_guess**(4./3.)))
            Sf_i2_new = ((mann(converge)**2)*guessQ(ii+1)*abs(guessQ(ii+1))/(A_i2_guess*A_i2_guess*Rh_i2_guess**(4./3.)))
            Sf_i_old = ((mann(converge)**2)*oldQ(ii)*abs(oldQ(ii))/(A_i_old*A_i_old*Rh_i_old**(4./3.)))
            Sf_i2_old = ((mann(converge)**2)*oldQ(ii+1)*abs(oldQ(ii+1))/(A_i2_old*A_i2_old*Rh_i2_old**(4./3.)))
            
            !Mo_i = guessQ(ii) + guessQ(ii+1) - oldQ(ii) - oldQ(ii+1) + &
            !    c*(theta(converge)*((guessQ(ii+1)*guessQ(ii+1)/A_i2_guess) - (guessQ(ii)*guessQ(ii)/A_i_guess)) + &
            !    (1-theta(converge))*((oldQ(ii+1)*oldQ(ii+1)/A_i2_old) - (oldQ(ii)*oldQ(ii)/A_i_old))) + &
            !    (2*delt*g)*(theta(converge)*0.5*(A_i_guess+A_i2_guess)*(((guessy(ii+1)-guessy(ii))/delx(converge)) -&
            !    So(converge) + 0.5*Sf_i_new) +&
            !    (1-theta(converge))*0.5*(A_i_old+A_i2_old)*(((oldy(ii+1)-oldy(ii))/delx(converge)) - So(converge) + 0.5*Sf_i_old))
            
            Mo_i = guessQ(ii) + guessQ(ii+1) - oldQ(ii) - oldQ(ii+1) + &
                c*theta(converge)*((guessQ(ii+1)*guessQ(ii+1)/A_i2_guess) - (guessQ(ii)*guessQ(ii)/A_i_guess) + &
                g*cen_i2_guess - g*cen_i_guess) + &
                c*(1-theta(converge))*((oldQ(ii+1)*oldQ(ii+1)/A_i2_old) - (oldQ(ii)*oldQ(ii)/A_i_old) + &
                g*cen_i2_old - g*cen_i_old) - &
                delt*theta(converge)*g*((So(converge) - Sf_i2_new)*A_i2_guess + (So(converge) - Sf_i_new)*A_i_guess) - &
                delt*(1-theta(converge))*g*((So(converge) - Sf_i2_old)*A_i2_old + (So(converge) - Sf_i_old)*A_i_old)
        
            ! Co. of  Momentum Derivatives
            Sf_yi = (mann(converge)**2)*guessQ(ii)*abs(guessQ(ii))*(((-4./3.)/((A_i_guess**2)*(Rh_i_guess**(7./3.))))*Rh_hi - &
                (2/((A_i_guess**3)*(Rh_i_guess**(4./3.))))*A_hi)
            Sf_yi2 = (mann(converge)**2)*guessQ(ii+1)*abs(guessQ(ii+1))*(((-4./3.)/((A_i2_guess**2)*(Rh_i2_guess**(7./3.))))*Rh_hi2 - &
                (2/((A_i2_guess**3)*(Rh_i2_guess**(4./3.))))*A_hi2)
            Sf_Qi = ((2*(mann(converge)**2)*(guessQ(ii)**2))/((A_i_guess**2)*(Rh_i_guess**(4./3.))*abs(guessQ(ii))))
            Sf_Qi2 = ((2*(mann(converge)**2)*(guessQ(ii+1)**2))/((A_i2_guess**2)*(Rh_i2_guess**(4./3.))*abs(guessQ(ii+1))))
            
            !Mo_yi = c*theta(converge)*((guessQ(ii)**2)/(A_i_guess**2))*A_hi + &
            !    (2*delt*g)*((theta(converge)/2)*(((guessy(ii+1)-guessy(ii))/delx(converge)) - So(converge) + 0.5*Sf_i_new)*A_hi + &
            !    theta(converge)*0.5*(A_i_guess+A_i2_guess)*((-1/delx(converge)) + 1.0*Sf_yi))
            !Mo_yi2 = -c*theta(converge)*((guessQ(ii+1)**2)/(A_i2_guess**2))*A_hi2 + &
            !    (2*delt*g)*((theta(converge)/2)*(((guessy(ii+1)-guessy(ii))/delx(converge)) - So(converge) + 0.5*Sf_i_new)*A_hi2 + &
            !    theta(converge)*0.5*(A_i_guess+A_i2_guess)*((1/delx(converge)) + 1.0*Sf_yi2))
            !Mo_Qi = 1 - (2*c*theta(converge)*guessQ(ii)/A_i_guess) - (2*delt*g)*(theta(converge)*0.5*(A_i_guess+A_i2_guess)*1.0*(Sf_Qi)) 
            !Mo_Qi2 = 1 + (2*c*theta(converge)*guessQ(ii+1)/A_i2_guess) - (2*delt*g)*(theta(converge)*0.5*(A_i_guess+A_i2_guess)*1.0*(Sf_Qi2)) 
            
            Mo_Qi = 1 - (2*c*theta(converge)*guessQ(ii)/A_i_guess) + delt*g*theta(converge)*A_i_guess*Sf_Qi
            Mo_Qi2 = 1 + (2*c*theta(converge)*guessQ(ii+1)/A_i2_guess) + delt*g*theta(converge)*A_i2_guess*Sf_Qi2
            Mo_yi = c*theta(converge)*((guessQ(ii)**2)/(A_i_guess**2))*A_hi - c*theta(converge)*g*cen_hi - &
                delt*theta(converge)*g*(A_hi*(So(converge) - Sf_i_new) + A_i_guess*(-Sf_yi))
            Mo_yi2 = -c*theta(converge)*((guessQ(ii+1)**2)/(A_i2_guess**2))*A_hi2 + c*theta(converge)*g*cen_hi2 - &
                delt*theta(converge)*g*(A_hi2*(So(converge) - Sf_i2_new) + A_i2_guess*(-Sf_yi2))

            if (Sf_Qi.ne.Sf_Qi) then  
                Mo_Qi = 1 - (2*c*theta(converge)*guessQ(ii)/A_i_guess)
            end if
            if (Sf_Qi2.ne.Sf_Qi2) then  
                Mo_Qi2 = 1 + (2*c*theta(converge)*guessQ(ii+1)/A_i2_guess)
            end if
        
            ! Edited Jacobian
            Je(2*ii,1) = Mass_Qi
            Je(2*ii,2) = Mass_yi
            Je(2*ii,3) = Mass_Qi2
            Je(2*ii,4) = Mass_yi2
            
            Je((2*ii)+1,1) = Mo_Qi
            Je((2*ii)+1,2) = Mo_yi
            Je((2*ii)+1,3) = Mo_Qi2
            Je((2*ii)+1,4) = Mo_yi2

            ! Residual Vector
            resid(2*ii) = -Mass_i
            resid((2*ii)+1) = -Mo_i
            
        end do
        
        deallocate(ql)
     
    end
