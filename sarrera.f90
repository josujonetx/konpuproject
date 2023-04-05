program sarrera

use tipoak

real(kind=dp), dimension(100):: x, r, v, a !100 partikula izango direlako 
integer, dimension(100):: c !kargen balioak
integer:: i, j,l
Real(kind=dp):: V0,U0,vm,w, dt, m

open(unit=111, file="K.dat", status="replace", action="write")

!sistemaren energia
U0=10.0

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

!----------------------------------------------------------------------------------------
!Txekeoa::

write(unit=*, fmt=*) sqrt((U0-V0)/50),"=", sqrt(sum(v*v))

!----------------------------------------------------------------------------------------

!x lista ordenatuko dugu errazagoa izan dadin ondorengo kalkuloa
do i=1,99 
    do j=i+1,100 
        if (x(i)>x(j)) then 
            m=x(i); x(i)=v(j); v(j)=m 
            if (j<50) then
                c(i)=1
            else (j>50) then
                c(i)=-1
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
do j=1,10000
    call karak(x,v,a,dt, l)
    do i=1,100
        x(i)=x(i)+v(i)*dt+a(i)/2*dt**2
        v(i)=v(i)+a(i)*dt
    enddo
    if (l=0) then
        v(1)=-v(1)
    else if (l=100) then
        v(100)=-V(100)
    else:
        m=x(l); x(l)=x(l+1); x(l+1)=m; m=v(l); v(l)=v(l+1); v(l+1)=m
    write(unit=111, fmt=*) dt*j, sum(v*v)/100/2
enddo


!__________________
contains
    function sinu(i,j) !potentziala + edo - izango den
    
    integer, intent(in)::i,j
    integer:: sinu
    
    sinu=(-1)**((i-1)/50)*(-1)**((j-1)/50)

    end function sinu
    
    subroutine karak(x,v,a,dt,l) !Denbora karakterestikoa
        
        real(kind=dp), dimension(:), intent(in):: x, a, v
        real(kind=dp), dimension(99)::m
        integer:: i,j
        real(kind=dp):: d,b,c,t1,t2
        real(kind=dp), intent(inout)::dt
        integer, intent(inout):: l
        
        ! 1. partikula ezker paretarekin talka 
        
        t1=(-x(1)-sqrt(v(1)**2-2*a(1)*x(1)))/a(1)
        t2=(-x(1)+sqrt(v(1)**2-2*a(1)*x(1)))/a(1)
        
        if ((t1<0.0_dp) .and. (t2>0.0_dp)) then
           dt=t2
           l=0
        else if (t1>0.0_dp) then
           dt=t1
           l=0
        end if
        
        !
        t1=(-x(1)-sqrt(v(1)**2-2*a(1)*x(1)))/a(1)
        t2=(-x(1)+sqrt(v(1)**2-2*a(1)*x(1)))/a(1)
        
        if ((t1<0.0_dp) .and. (t2>0.0_dp)) then
           dt=t2
           l=0
        else if (t1>0.0_dp) then
           dt=t1
           l=0
        end if
        
        t1=(1-x(1)-sqrt(v(1)**2-2*a(1)*(1-x(1))))/a(1)
        t2=(1-x(1)+sqrt(v(1)**2-2*a(1)*(1-x(1))))/a(1)
        if ((t1<0.0_dp) .and. (t2>0.0_dp)) then
           karak=t2
           l=100
        else if (t1>0.0_dp) then
           karak=t1
           l=100
        end if
        
        do i=1,100
            c=x(i+1)-x(i)
            b=v(i+1)-v(i)
            d=(a(i+1)-a(i))/2
            t1=(-b+sqrt(b**2-4*d*c))/2/d
            t2=(-b+sqrt(b**2-4*d*c))/2/d
        if ((t1<0.0_dp) .and. (t2>0.0_dp)) then
           dt=t2
           l=i
        else if (t1>0.0_dp) then
           dt=t1
           l=i
        end if
        enddo
    end subroutine karak
    
end program sarrera
