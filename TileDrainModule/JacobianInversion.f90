subroutine JacobianInversion()

    use global
    
    implicit none

    double precision, dimension(2*N(converge)) :: m3, m4, m2, z
    
    m3(1) = Je(1,3)
    m4(1) = Je(1,4)
    z(1) = resid(1)

    do i = 2,2*N(converge)
        if (mod(i,2).eq.0) then
            m2(i) = -(Je(i,1)*m4(i-1)/m3(i-1)) + Je(i,2)
            z(i) = -(Je(i,1)*z(i-1)/m3(i-1)) + resid(i)
        else
            m2(i) = -(Je(i,1)*m4(i-2)/m3(i-2)) + Je(i,2)
            m3(i) = -(Je(i-1,3)*m2(i)/m2(i-1)) + Je(i,3)
            m4(i) = -(Je(i-1,4)*m2(i)/m2(i-1)) + Je(i,4)
            z(i) = -(m2(i)*z(i-1)/m2(i-1)) - (Je(i,1)*z(i-2)/m3(i-2)) + resid(i)
        end if 
    end do

    
    x(2*N(converge)) = z(2*N(converge))/m2(2*N(converge))    

    do i = (2*N(converge)-1),1,-1
        if (mod(i,2).ne.0) then
            x(i) = (z(i)-(m4(i)*x(i+1)))/m3(i)
        else 
            x(i) = (z(i) - Je(i,4)*x(i+2) - (Je(i,3)*x(i+1)))/m2(i)
        end if
    end do


    end