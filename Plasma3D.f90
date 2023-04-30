module Plasma3D

public:: Oreka


contains
    function Oreka(n) 
        use tipoak
        use funtzioak
        use rk4
        
        real(kind=dp), intent(in):: n
        integer             :: i,j, l,o
        real(kind=dp)       :: t, h, U0,V0,vm,w
        real(kind=dp), dimension(n)  :: x,y,z,c,r, theta, phy,vx,vy,vz
        real(kind=dp), dimension(6*n) :: p
        real(kind=dp), dimension(1000,6*n):: Oreka
        integer, parameter  :: m=1000
        real(kind=dp), parameter :: ta=0.0_dp, tb=100.0_dp


        h=(tb-ta)/m

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

        do i=1,n
        do j=1,n
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

        do i=1,n
            p(6*i-5)= x(i)
            p(6*i-4)= y(i)
            p(6*i-3)= z(i)
            p(6*i-2)= vx(i)
            p(6*i-1)= vy(i)
            p(6*i)  = vz(i)
        enddo

        do i=1,m
                call rk4_paso_dp(t,p, Plasma ,h)
                (Oreka(i,j)=p(j), j=1,6*n)
        end do

        end function Oreka

end module Plasma3D
