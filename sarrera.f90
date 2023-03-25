program plasma
use tipoak

real(kind=dp), dimension(100):: x, v !100 partikula izango direlako 
integer:: i, j 
Real(kind=dp):: V0

call random_number(x) !100 posizio aleaorio


!Potentzialaren kalkuloa:

V0=0.0_dp
do i=1,100
do j=1,100
    V0=V0+sinu(i,j)/abs(x(i)-x(j))
enddo
enddo
V0=V0/2

!__________________
contains
    function sinu(i,j) !potentziala + edo - izango den
    integer, intent(in)::i,j
    integer:: sinu
    sinu=(-1)**(floor((i-1)/50)*floor((j-1)/50))
    end function sinu
end program plasma
