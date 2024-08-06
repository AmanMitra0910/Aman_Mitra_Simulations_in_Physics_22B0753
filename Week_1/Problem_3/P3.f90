program P3
    implicit none
    real :: t_end , h1 , h2 

    t_end = 10.0
    h1 = 0.1
    h2 = 0.01

    call euler(h1, t_end)
    call euler(h2 , t_end)

contains

    subroutine euler(h , t_end)
        implicit none
        real, intent(in) :: h , t_end
        integer :: n_steps , i
        real :: t , y , v , y1 , v1 , y_exact
        character(len = 100) :: estimated_euler
        character(len = 10) :: h_str

        write(h_str, '(f5.3)') h  
        estimated_euler = 'output_h_' // trim(h_str) // '.txt'
        n_steps = int(t_end / h)

        open(unit = 10, file =estimated_euler, status='replace')

        y = 0.0
        v = 5.0
        t = 0.0

        write(10, *) "   t   " , "  y_estimated_euler  " , "  y_exact  "

        do i = 0, n_steps
            y_exact = exact(t)

            write(10 ,*) t , y , y_exact

            y1 = y + h*v
            v1 = v + h*(10.0 * sin(t) - 5.0 * v - 6.0 * y)

            y = y1
            v = v1
            t = t + h
        
        end do

        close(10)

        print *, "Results for h = ", h ," are documented in" , trim(estimated_euler)
    
    end subroutine euler

    function exact(t) result(y_exact)
        implicit none
        real, intent(in) :: t
        real :: y_exact

        y_exact = -6.0 *exp(-3.0 * t) + 7.0 * exp(-2.0 * t) + sin(t) - cos(t)

    end function exact

end program P3