program Q43
    implicit none
    real :: x, sum
    integer :: n

    print *, "Enter a number:"
    read(*,*) x
    print *, "Enter the number of iterations:"
    read(*,*) n
    sum = series_sum(x , n)

    print *, "The sum of the series is", sum

    print *, " To verify the value of cos(x) using this operation  is, if x = 1 ," , series_sum(1.0 , 1000000)
    print *, " if x = 4 ,", series_sum(4.0 , 1000000)

contains

    real function series_sum(x , n)
        implicit none
        integer , intent(in) :: n
        real , intent(in) :: x
        integer :: i
        real :: term, sum

        sum = 0.0

        do i = 0, n-1
            term = ((-1)**i)*(x**(2*i))/fac(2*i)
            sum = sum + term
        end do

        series_sum = sum
    end function series_sum
    
    real function fac(n)
        implicit none
        integer , intent(in) :: n
        integer :: i
        real :: a

        if ( n == 0) then
            fac = 1
        else    
            i = n
            a = 1
            do while (i>=1)
               a = a*i
               i = i - 1
            end do 

            fac = a
        end if

    end function fac
        
end program Q43