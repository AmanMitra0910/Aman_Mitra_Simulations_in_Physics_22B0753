program Q2
   implicit none
   
   real :: input(3)
   real :: a , func
   integer :: i

   input = (/ -0.5 , 0.1 , 2.5 /)

   do i = 1, 3
      a = input(i)
      if (a<2.0) then
        func = 5.0*a*a + 3.0*a + 2.0
      else if (a==2.0) then
        func = 0.0
      else if (a>2.0) then
        func = 5.0*a*a - 3.0*a + 1.0
      end if   

      print *, "if input is" , a ,"then, the value of the function is" , func
   end do

end program Q2

      
