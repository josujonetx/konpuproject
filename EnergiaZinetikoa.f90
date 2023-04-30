program EnergiaZinetikoa

use Plasma3D

real(kind=dp), dimension(100,1000):: K,V0
integer, parameter:: n=50 !Partikula kop.
real(kind=dp), dimension(1000, 6*n)::p
real(kind=dp), dimension(n):: x,y,z,vx,vy,vz
real(kind=dp), dimension(1000):: km,vm
integer:: i,j,l,o

open(unit=11, file="EnergiaZinetikoa.dat", status="replace", action="write")

do i=1,100
   p=Oreka(n)
   do l=1,1000
      do j=1,n
         x(j)=p(l,6*j-5)
         y(j)=p(l,6*j-4)
         z(j)=p(l,6*j-3)
         vx(j)=p(l,6*j-2)
         vy(j)=p(l,6*j-1)
         vz(j)=p(l,6*j)
      enddo
      K(i,l)=sum(vx**2+vy**2+vz**2)/2
      V0(i,l)0.0_dp
      do j=1,n
        do o=1,n
            if (o/=j) then
              V0(i,l)=V0(i,l)-c(o)*c(j)/sqrt((x(o)-x(j))**2+(y(o)-y(j))**2+(z(o)-z(j))**2)
            end if
         enddo
      enddo
   enddo
enddo


do i=1,100
   km=km+K(i,:)/100
   vm=vm+V(i,:)/100
enddo


do l=1,1000
   write(unit=11, fmt=*) km(l), vm(l), km(l)+vm(l)
enddo

close(11)

end program EnergiaZinetikoa
