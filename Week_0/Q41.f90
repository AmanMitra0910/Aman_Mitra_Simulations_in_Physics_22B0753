program Q41
    implicit none
    real :: x, sum
    integer :: n

    print *, "Enter a number:"
    read(*,*) x
    print *, "Enter the number of iterations:"
    read(*,*) n
    sum = series_sum(n)

    print *, "The sum of the series is", sum

    print *, " To verify the value of pi(3.141592653) using this operation  is," , series_sum(1000000000)

contains

    real function series_sum(n)
        implicit none
        integer , intent(in) :: n
        integer :: i
        real :: term, sum

        sum = 0.0

        do i = 0, n-1
            term = 4.0 / (2*i + 1) * (-1) ** i
            sum = sum + term
        end do

        series_sum = sum
    end function series_sum

end program Q41
