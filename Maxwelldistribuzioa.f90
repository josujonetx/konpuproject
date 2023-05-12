program maxwell

use tipoak
use sarrera

integer, parameter::q=60, m=100 !m=Errepikapen kop., n=Partikula kop., o=Iterazio kop.
real(kind=dp), dimension(o, n+1 ,4)::p
real(kind=dp), dimension(m*n):: v
integer:: i,j,l
integer, dimension(q):: k
real(kind=dp):: v1max, v1min, vmax,vmin, a, d, dv, vm

open(unit=111, file="Abiad_distr.dat", status="replace",action="write" )
open(unit=112, file="Maxwell.dat", status="replace",action="write" )

do j=0,m-1
    call Plasma1D(p)
    do i=1,n
       v(n*j+i)=p(o,i+1,2)
    enddo

enddo

    do i=1,n*m-1
      do l=i+1,n*m
        if (v(i)>v(l)) then
            a=v(i)
            v(i)=v(l)
            v(l)=a
         end if
       end do
     end do
    v1max=v(n*m)
    v1min=v(1)
    d=v1max-v1min
    do i=1,q
      k(i)=0
      vmax=v1min+d/q*i
      vmin=v1min+d/q*(i-1)
      do l=1,m*n
        if ((vmin<=v(l)) .and. (v(l)<=vmax)) then
          k(i)=k(i)+1
        end if
      enddo
    enddo
    write(unit=*, fmt=*) sum(k)

do i=1,q
        write(unit=111, fmt=*) v1min+d/q*(i-1), real(k(i))/n*q/d/m
        write(unit=111, fmt=*) v1min+d/q*(i), real(k(i))/n*q/d/m
enddo
vm=sum(v*v)/size(v)
do i=1,1000
        dv=v1min+d/1000*(i-0.5)
        write(unit=112, fmt=*) dv, 1.0_dp/sqrt(2*acos(-1.0)*vm)*Exp(-dv**2/2/vm)
enddo

end program maxwell
