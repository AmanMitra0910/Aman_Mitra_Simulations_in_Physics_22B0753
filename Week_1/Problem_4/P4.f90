program P4
    implicit none
    real :: h
    
    h = 0.1

    call RK4(h)

contains
    
    subroutine RK4(h)
        implicit none
        real, intent(in) :: h
        integer :: i , n
        real :: t , y , v , y1 , v1 , y_exact
        character(len = 100) :: estimated_RK4
        character(len = 10) :: h_str

        write(h_str, '(f5.3)') h  
        estimated_RK4 = 'output_h_' // trim(h_str) // '.txt'

        open(unit = 10, file =estimated_RK4, status='replace')

        n = 200
        t = 0

        write(10, *) "   t   " , "  y_estimated_RK4  " , "  y_exact  "

        do i = 0, n
            y_exact = exact(t)

            write(10 ,*) t , y , y_exact

            y1 = f_y(t , y , v)
            v1 = f_v(t , y , v)
            
            y = y1
            v = v1
            t = t + h/2
        
        end do

        close(10)

        print *, "Results for h = ", h ," are documented in" , trim(estimated_RK4)

    end subroutine RK4

    function f_y(t , y , v) result(y1)
        implicit none
        real, intent(in) :: t , y , v
        integer :: i
        real :: k(4) , l(4) , y1 , h , lv , lt , ly
        
        h = 0.1
        lt = t
        lv = v
        ly = y

        do i = 1 , 3
            
            k(i) = lv
            l(i) = 10 * sin(t) - 6 * ly - 5 * lv

            ly = ly + k(i) * (h/2)
            lv = lv + l(i) * (h/2)
            lt = lt + h/2
        
        end do
        
        k(4) = ly + l(3) * (h/2)

        y1 = y + (h/6) * (k(1) + 2*k(2) + 2*k(3) + k(4))
    
    end function f_y

    function f_v(t , y , v) result(v1)
        implicit none
        real, intent(in) :: t , y , v
        integer :: i
        real :: k(4) , l(4) , v1 , h , lv , lt , ly

        h = 0.1
        lt = t
        lv = v
        ly = y

        do i = 1 , 3
            
            k(i) = lv
            l(i) = 10 * sin(t) - 6 * ly - 5 * lv

            ly = ly + k(i) * (h/2)
            lv = lv + l(i) * (h/2)
            lt = lt + h/2
        
        end do

        l(4) = 10 * sin(lt + (h/2)) - 6 * (ly + k(3) * (h/2)) - 5 * (lv + l(3) * (h/2))

        v1 = v + (h/6) * (l(1) + 2*l(2) + 2*l(3) + l(4))
    
    end function f_v

    function exact(t) result(y_exact)
        implicit none
        real, intent(in) :: t
        real :: y_exact

        y_exact = -6.0 *exp(-3.0 * t) + 7.0 * exp(-2.0 * t) + sin(t) - cos(t)

    end function exact

end program P4    