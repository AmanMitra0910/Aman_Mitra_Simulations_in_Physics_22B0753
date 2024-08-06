program Q12
   implicit none

   real :: c
   real :: d
   real :: a
   real :: m
   real :: s
   real :: result
   real :: ten 
   
   ten =10
   c = 1.0
   d = 2.0
   a = 15
   m = (c**c)+(d**d)
   s = c*d
   result = (log(c)/log(ten)) + cosd(a) + abs(m) + 2*sqrt(s)
   
   print *, "the answer to the operation is" , result

end program Q12   