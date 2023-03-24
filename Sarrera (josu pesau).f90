program plasma
use tipoak

real(kind=dp), dimension(100):: x, v !100 partikula izango direlako 
integer:: i, j 
Real(kind=dp):: V0
call random_number(x)
call random_number(v)

do i=1,100
do j=1,100
    V0=V0+((50-i)*(50-j))/abs((50-i)*(50-j))/abs(x(i)-x(j))
enddo
enddo
end program plasma