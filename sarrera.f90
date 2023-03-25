program sarrera

use tipoak

real(kind=dp), dimension(100):: x, r, v, E !100 partikula izango direlako 
integer:: i, j,l
Real(kind=dp):: V0,U0,vm,w, dt,k

open(unit=111, file="K(t).dat", status="replace", action="write")

!sistemaren energia
U0=10.0

potentziala: do
!Potentzialaren kalkuloa:

call random_number(x) !100 posizio aleaorio

V0=0.0_dp

do i=1,100
do j=1,100
    V0=V0-sinu(i,j)*abs(x(i)-x(j))
enddo
enddo

!kalkuluak simplifikatzeko interakzio guztiak 2 bider batu dira


V0=V0/2 

!potentzialak ezin duenez sistemaren energia gainditu prozesua errepikatuko da:

if (V0<U0) then
    exit potentziala
end if

end do potentziala
! Partiukula ez erlatibistak (E=1/2·m·v^2) eta eta denek m=1 onartuko dugu. 
! Hortaz abiaduraren bataz besteko kuadratikoa-ren erroa:

vm=sqrt((U0-V0)/50)

call random_number(r) !abiaduren aleatoriotasuna bermatzeko

w=sqrt(sum(r*r))

v=r*vm/w !100 elementuko abiaduren lista bat

!----------------------------------------------------------------------------------------
!Txekeoa::

write(unit=*, fmt=*) sqrt((U0-V0)/50),"=", sqrt(sum(v*v))

!----------------------------------------------------------------------------------------

!eboluzio temporala 0.01-eko jauziez kalkulatua
dt=0.01 
do l=1,600
!Eremu elektrikoa::
do j=1,100
    E(j)=0.0_dp
    do i=1,100
        E(j)=E(j)+sig(i,j)
    enddo
enddo
v=v+E*dt
x=x+v*dt+E/2*dt**2

K=sqrt(sum(v*v))/100/2

enddo
!__________________
contains
    function sinu(i,j) !potentziala + edo - izango den
    
    integer, intent(in)::i,j
    integer:: sinu
    
    sinu=(-1)**((i-1)/50)*(-1)**((j-1)/50)

    end function sinu
    
    function sig(i,j) !eskubi ala ezker
    
    integer, intent(in)::i,j
    integer:: sig
    
    if (i==j) then
        sig=0
    else if ((x(j)-x(i))>0) then
        sig=1
    else: 
        sig=-1
    end if
    
    end function sig
    
end program sarrera
