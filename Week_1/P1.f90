program P1
    implicit none
    real(8) :: x , L , m , h , E , a , k , ee , fx , dfx
    integer :: i

    x = 1.0
    L = 2.0e-10
    m = 9.1e-31
    h = 6.58e-16
    ee = 1.6e-19
    a = 1.0e-10

    do i = 1 , 10000
        
        call deriv(x , fx , dfx)

        x = x - (fx / dfx)

        if (abs(fx) < a) exit
    
    end do

    k = x/L
    E = (k*h)**2 / (2*m*ee)

    print *, "Calculated Ground state energy is ,", E ,"eV"

contains

    subroutine deriv(x, fx, dfx)
        real(8), intent(in) :: x
        real(8), intent(out) :: fx , dfx
        fx = cos(x) - 0.226*x
        dfx = -sin(x) - 0.226
    end subroutine deriv   
    
end program P1    
