program Q3
    implicit none
    integer :: N , i , j , a , b
   
    a = 1
    b = 1
    print *, "Enter  a number: "
    read(*,*) N
    if (N == 0) then
      print *, "the factorial of the following number is, 1"
    else    
      i = N
      do while (i>=1)
         a = a*i
         i = i - 1
      end do 

      print *, "the factorial of the following number is," , a
    end if
    
    j = N
    if (mod(N , 2) == 0) then
      print *, "the double factorial of this number is not possible."
    else 
      do while (j>=1)
        b = b*j
        j = j - 2   
      end do
      print *, "the double factorial of this number is," , b
    end if      
      
end program Q3    