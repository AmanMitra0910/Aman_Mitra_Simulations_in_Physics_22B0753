program Q11
  implicit none

  real :: a
  real :: b
  real :: exp
  real :: sub
  real :: result
  
  a = 3.0
  b = 4.0
  exp = a**b
  sub = b-a
  result = exp/(b*sub)

  print *, "the answer to this operation is", result

end program Q11    