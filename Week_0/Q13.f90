program Q13
   implicit none

   real :: e
   real :: pi
   real :: c
   real :: m
   real :: a
   real :: result

   e = 2.71828
   pi = 3.1415926
   c = 2.0
   m = 1.0
   a = 2.0
   result =(e**(sqrt(2*a*((c-m)**3))))/(a*(sqrt(2*pi)))

   print *, "the answer to this operation is", result

end program Q13   