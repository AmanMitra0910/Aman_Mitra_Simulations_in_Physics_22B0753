program P1
    implicit none
    real :: d(3 , 3) , lu(3 , 3) , x(1 , 3) , b (1 , 3) , luxx(1 , 3) , s(1 , 3) , di(3 , 3) , y(1 , 3) , a(3 , 3)
    real :: t
    integer :: i , n_max

    d = reshape([4 ,0 ,0 ,0 ,6 ,0 ,0 ,0 ,7], shape(d))
    lu = reshape([0 ,-1 ,-1 ,-2 ,0 ,1 , -1, 1, 0], shape(lu))
    b = reshape([3 ,9 ,-6], shape(b))
    x = reshape([0, 0, 0], shape(x))
    a = 1e-5
    n_max = 100

    do i = 1 , n_max
        
        luxx = mlt13(lu , x)
        s = b - luxx
        di = inv33(d)
        a = d + lu

        x = mlt13(di , s)
        y = mlt13(a , x) - b

        if((y(1 , 1) < t) .AND. (y(1 , 2) < t) .AND. (y(1 , 3) < t)) exit
    end do

    print *, "The solution for this set of linear equations using Jacobi Method is, x =", x(1 , 1), " y =", x(1, 2), " z =", x(1, 3)

contains

    real function mlt13(m1 , m2)
        implicit none
        real, intent(in) :: m1(3 , 3) , m2(1 , 3)
        real :: mlt13(1 , 3)
        integer :: i

        do i = 1, 3
            mlt13(1, i) = sum(A(i, :) * B(1, :))
        end do
    end function mlt13

    real function inv33(m) 
        real, intent(in) :: m(3 , 3)
        real :: inv33(3 , 3)
        real :: adj(3 , 3)
        
        adj(1,1) =  m(2,2)*m(3,3) - m(2,3)*m(3,2)
        adj(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
        adj(1,3) =  m(1,2)*m(2,3) - m(1,3)*m(2,2)
        adj(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
        adj(2,2) =  m(1,1)*m(3,3) - m(1,3)*m(3,1)
        adj(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
        adj(3,1) =  m(2,1)*m(3,2) - m(2,2)*m(3,1)
        adj(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
        adj(3,3) =  m(1,1)*m(2,2) - m(1,2)*m(2,1)

        inv33 = transpose(adj) / det(m)

    end function inv33

    real function det(m)
        real, dimension(3, 3) :: m
        det = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))

    end function det

end program P1