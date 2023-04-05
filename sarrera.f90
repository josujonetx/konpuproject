program sarrera

use tipoak

real(kind=dp), dimension(100):: x, r, v, a !100 partikula izango direlako 
integer, dimension(100):: c !kargen balioak
integer:: i, j,l,k
Real(kind=dp):: V0,U0,vm,w, dt, m,t

open(unit=111, file="K.dat", status="replace", action="write")

!sistemaren energia
U0=10.0_dp

potentziala: do
!Potentzialaren kalkuloa:

call random_number(x) !100 posizio aleaorio [0,1)

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

!------------------------------------------------------------------------------>
!Txekeoa::

write(unit=*, fmt=*) sqrt((U0-V0)/50),"=", sqrt(sum(v*v))

!------------------------------------------------------------------------------>

!x lista ordenatuko dugu errazagoa izan dadin ondorengo kalkuloa
do i=1,99
    do j=i+1,100
        if (x(i)>x(j)) then
            m=x(i)
            x(i)=v(j)
            v(j)=m
            if (j<=50) then
                c(i)=1   !kargen balioa ordena aldatzean ez baitauzkagu hasierako 50na + eta  - 
            else
                c(i)=-1
            end if
        end if
    end do
end do

! Azelerazioa:
do i=1,100
    a(i)=0.0_dp
    do j=1,i-1
        a(i)=a(i)+c(i)*c(j)
    enddo
    do j=i+1,100
        a(i)=a(i)-c(i)*c(j)
    enddo
enddo

!dinamika

do j=1,10000
    call karak(x,v,a,dt, l) !Denbora karakteristikoa (dt) eta prozesua (l)
    do i=1,100
        x(i)=x(i)+v(i)*dt+a(i)/2*dt**2   ! Eboluzio denborala
        v(i)=v(i)+a(i)*dt
    enddo
    if (l==0) then
        v(1)=-v(1) !ezkerrako paretaren aurkako talka
    else if (l==100) then
        v(100)=-V(100) !eskuineko paretaren aurkako talka
    else !gurutzatzea
        m=x(l)
        x(l)=x(l+1)
        x(l+1)=m
        m=v(l)
        v(l)=v(l+1)
        v(l+1)=m
        m=a(i)
        a(i)=a(i)+2*c(l)*c(l+1)
        a(i+1)=a(i+1)-2*c(l)*c(l+1)
        m=c(i)
        c(i)=c(i+1)
        c(i+1)=m
    end if
    t=t+dt
    write(unit=111, fmt=*) t, sum(v*v)/100/2
    
    !txekeoa??----------------------------------------------------    
    do i=1,100
        write(unit=*, fmt=*) "Xmax", X(100)
    enddo
    V0=0.0_dp
    do i=1,100
    do k=1,100
        V0=V0-c(i)*c(k)*abs(x(i)-x(k))
    enddo
    enddo
    V0=V0/2
    write(unit=*, fmt=*) sqrt((U0-V0)/50),"=", sqrt(sum(v*v))
    
    !----------------------------------------------------
    
enddo

contains
    function sinu(i,j) !potentziala + edo - izango den

    integer, intent(in)::i,j
    integer:: sinu

    sinu=(-1)**((i-1)/50)*(-1)**((j-1)/50)

    end function sinu

    subroutine karak(x,v,a,dt,l) !Denbora karakterestikoa

        real(kind=dp), dimension(:), intent(in):: x, a, v
        integer:: i,j
        real(kind=dp):: d,b,c,t1,t2
        real(kind=dp), intent(inout)::dt
        integer, intent(inout):: l
        
            ! 1. partikula ezker paretarekin talka 

        t1=(-v(1)-sqrt(v(1)**2-2*a(1)*x(1)))/a(1)
        t2=(-v(1)+sqrt(v(1)**2-2*a(1)*x(1)))/a(1)

        if ((t1<0.0_dp) .and. (t2>0.0_dp)) then
           dt=t2
           l=0
        else if (t1>0.0_dp) then
           dt=t1
           l=0
        end if
            
             ! 100. partikula eskuin paretarekin talka 
             
        t1=(-v(100)-sqrt(v(100)**2-2*a(100)*(-1.0+x(100))))/a(100)
        t2=(-v(100)+sqrt(v(100)**2-2*a(100)*(-1.0+x(100))))/a(100)
        if ((t1<0.0_dp) .and. (t2>0.0_dp) .and. (dt>t2)) then
           dt=t2
           l=100
        else if (t1>0.0_dp .and. (dt>t1) ) then
           dt=t1
           l=100
        end if
            
            !Partikulak elkar gurutzatu
            
        do i=1,100
            c=x(i+1)-x(i)
            b=v(i+1)-v(i)
            d=(a(i+1)-a(i))/2
            t1=(-b-sqrt(b**2-4*d*c))/2/d
            t2=(-b+sqrt(b**2-4*d*c))/2/d
            if ((t1<0.0_dp) .and. (t2>0.0_dp) .and. (dt>t2) ) then
               dt=t2
               l=i
            else if (t1>0.0_dp .and. (dt>t1)) then
               dt=t1
               l=i
            end if
        enddo
    end subroutine karak
end program sarrera


