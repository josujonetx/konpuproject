module Plasma3D
use tipoak
use funtzioak
use rk4

integer             :: i,j, l,o
real(kind=dp)       :: t, h, U0,V0,vm,w
real(kind=dp), dimension(50)  :: x,y,z,c,r, theta, phy,vx,vy,vz
real(kind=dp), dimension(300) :: p
real(kind=dp), dimension(100,1000):: K
real(kind=dp), dimension(1000):: km
integer, parameter  :: m=1000
real(kind=dp), parameter :: ta=0.0_dp, tb=100.0_dp

open(unit=12,file="plasma3D",status="replace",action="write")

do j=1,100
    (k(j,i)=0.0_dp, i=1,1000)
enddo


h=(tb-ta)/m

do l=1,100
t = 0.0_dp


!sistemaren energia
U0=10.0_dp

potentziala: do
!Potentzialaren kalkuloa:


call random_seed()
call random_number(x,y,z) !50 posizio aleaorio [0,1)

do i=1,n/2 !kargen balioak 1/-1
    c(i)=1
    c(i+n/2)=-1
enddo

c(n)=-1 !partikula kopurua bakoitia balitz

V0=0.0_dp

do i=1,50
do j=1,50
    if (i/=j) then
      V0=V0-c(i)*c(j)/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
    end if
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

vm=sqrt((U0-V0)/n*2)

call random_number(r,theta,phy) !abiaduren aleatoriotasuna bermatzeko r-[0,1)

r=2*r-1 ! r->(-1,1)
theta=theta*acos(-1.0)*2 ! theta->[0,2pi)
phy=phy*acos(-1.0)       ! phy-> [0,pi)

w=sqrt(sum(r*r))

v=r*vm/w !50 elementuko abiaduren lista bat

vx=v*sin(phy)*cos(theta)
vy=v*sin(phy)*sin(theta)
vz=v*cos(phy)

do i=1,50
    p(6*i-5)= x(i)
    p(6*i-4)= y(i)
    p(6*i-3)= z(i)
    p(6*i-2)= vx(i)
    p(6*i-1)= vy(i)
    p(6*i)  = vz(i)
enddo

do i=1,m
        call rk4_paso_dp(t,p, Plasma ,h)
        V0=0.0_dp

        do j=1,50
        do o=1,50
            if (i/=o) then
              V0=V0-c(j)*c(o)/sqrt((x(o)-x(j))**2+(y(o)-y(j))**2+(z(o)-z(j))**2)
            end if
        enddo
        enddo
        V0=V0/2
        k(l,i)=sum( p(6*i-2)**2+ p(6*i-1)**2+ p(6*i)**2)/2
end do



enddo
close(unit=12)


end module Plasma3D
