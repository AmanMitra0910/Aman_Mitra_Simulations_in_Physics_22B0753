program P5
    implicit none
    real, parameter :: s = 10.0, r = 28.0, b = 8.0/3.0, t_end = 50.0, dt = 0.01
    integer, parameter :: steps = int(t_end / dt)

    real :: t
    real, dimension(3) :: y1 , y2 , dy_dt
    integer :: i

    y1 = (/1.0 , 1.0 , 1.0/)
    y2 = (/ 1.01 , 1.0 , 1.0/)

    open(unit = 1 , file = 'lorenz_sol1.txt')
    open(unit = 2 , file = 'lorenz_sol2.txt')

    do i = 1, steps
        t = (i-1) * dt

        call lorenz(t , y1 , dy_dt)
        y1 = y1 + dy_dt * dt
        write(1 , *) t , y1

        call lorenz(t , y2 , dy_dt)
        y2 = y2 + dy_dt * dt
        write(2 , *) t , y2

    end do  

    close(1)
    close(2)
    
contains

    subroutine lorenz(t , y , dy_dt)
        implicit none
        real, intent(in) :: t
        real, intent(in), dimension(3) :: y
        real, intent(out), dimension(3) :: dy_dt

        dy_dt(1) = s * (y(2) - y(1))
        dy_dt(2) = (y(1) * (r - y(3))) - y(2)
        dy_dt(3) = y(1) * y(2) - b * y(3)

    end subroutine lorenz

end program P5