program Q42
    implicit none
    real :: x, sum
    integer :: n

    print *, "Enter a number:"
    read(*,*) x
    print *, "Enter the number of iterations:"
    read(*,*) n
    sum = series_sum(x , n)

    print *, "The sum of the series is", sum

    print *, " To verify the value of e^x using this operation  is, if x = 1 ," , series_sum(1.0 , 1000000000)
    print *, " if x = 4 ,", series_sum(4.0 , 1000000000)

contains

    real function series_sum(x , n)
        implicit none
        integer , intent(in) :: n
        real , intent(in) :: x
        integer :: i
        real :: term, sum

        sum = 0.0

        do i = 0, n-1
            term = (x**i)/fac(i+1)
            sum = sum + term
        end do

        series_sum = sum
    end function series_sum
    
    real function fac(n)
        implicit none
        integer , intent(in) :: n
        integer :: i
        real :: a

        i = n
        a = 1
        do while (i>=1)
            a = a*i
            i = i - 1
        end do 

        fac = a
    end function fac
        
end program Q42