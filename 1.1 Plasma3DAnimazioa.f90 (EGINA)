program Animazioa

use tipoak
use Plasma3D
use funtzioak

real(kind=dp), dimension(m):: K,V0
integer, parameter:: n=100 !Partikula kop.
real(kind=dp), dimension(m, 6*n)::p
real(kind=dp), dimension(n):: x,y,z,vx,vy,vz,c 
integer:: i,j,l,o
real(kind=dp)::r

open(unit=111, file="Donuta.dat", status="replace", action="write")


 do i=1,n/2 !kargen balioak 1/-1
            c(i)=1.0_dp
            c(i+n/2)=-1.0_dp
 enddo

 c(n)=-1.0_dp !partikula kopurua bakoitia balitz


 p=Oreka(n)
 do l=1,m
     do j=1,n
         x(j)=p(l,6*j-5)
         y(j)=p(l,6*j-4)
         z(j)=p(l,6*j-3)
         vx(j)=p(l,6*j-2)
         vy(j)=p(l,6*j-1)
         vz(j)=p(l,6*j)
     enddo
     K(l)=sum(vx*vx+vy*vy+vz*vz)
     v0(l)=0.0_dp
     do j=1,n
        do o=1,n
            r=sqrt((x(o)-x(j))**2+(y(o)-y(j))**2+(z(o)-z(j))**2)
            if (o/=j) then
              V0(l)=V0(l)-c(o)*c(j)/r
            end if
         enddo
      enddo
enddo

do i=1,n
        write(unit=111, fmt=*) (p(l,6*i-5),p(l,6*i-4), p(l,6*i-3), l=1,m)
enddo

close(111)
end program Animazioa
