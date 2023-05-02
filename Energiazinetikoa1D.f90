program EnergiaZinetikoa

use tipoak
use Plasma3D

real(kind=dp), dimension(100,m):: K,V0
integer, parameter:: n=50 !Partikula kop.
real(kind=dp), dimension(m, 2*n)::p
real(kind=dp), dimension(n):: x,v
real(kind=dp), dimension(m):: km,vm
integer:: i,j,l,o

open(unit=11, file="EnergiaZinetikoa1D.dat", status="replace", action="write")

 do i=1,n/2 !kargen balioak 1/-1
            c(i)=1
            c(i+n/2)=-1
 enddo

 c(n)=-1 !partikula kopurua bakoitia balitz



do i=1,100
   p=Plasma1D(n)
   do l=1,m
      do j=1,n
         x(j)=p(l,2*j-1)
         v(j)=p(l,2*j)
      enddo
      K(i,l)=sum(v*v)/2
      V0(i,l)=0.0_dp
      do j=1,n
        do o=1,n
            if (o/=j) then
              V0(i,l)=V0(i,l)-c(o)*c(j)/((x(o)-x(j))**2
            end if
         enddo
      enddo
   enddo
enddo

do j=1,m
     km(j)=0.0_dp
     vm(j)=0.0_dp
enddo

do i=1,100
do j=1,m
   km(j)=km(j)+K(i,j)/100
   vm(j)=vm(j)+V0(i,j)/100
enddo
enddo

do l=1,m
   write(unit=11, fmt=*) l,km(l), vm(l), km(l)+vm(l)
enddo

close(11)

end program EnergiaZinetikoa
