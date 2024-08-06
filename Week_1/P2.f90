program P2 
    implicit none
    real :: x0 , x1 , fx0 , fx1 , tolerance , temp

    tolerance = 1.0e-20
    x0 = 1.0
    x1 = 3.0
    
    call function(x0 , fx0)
    call function(x1 , fx1)

    do while (abs(fx1) > tolerance)
        
        temp = x1
        x1 = x1 - ((x0 - x1) / (fx0 - fx1)) * fx1
        x0 = temp

        call function(x0 , fx0)
        call function(x1 , fx1)

    end do 
    
    print *, "The Root of the function is ," , x1 
    print *, "and the value of the function is ,", fx1

contains 

    subroutine function(x , fx)
        real, intent(in) :: X
        real, intent(out) :: fx
        fx = exp(2.0 * x) + x - 5.0
    end subroutine function

end program P2
