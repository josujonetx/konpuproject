program pant

use egoerak

real(kind=dp),dimension(49)::Pp,Pm
real(kind=dp),dimension(100,4)::m
real(kind=dp),dimension(99,2)::karg
integer::k,j,i,n
real(kind=dp)::d,q,p,s,l

open(unit=11,status="replace",action="write",file="Pp.dat")
open(unit=12,status="replace",action="write",file="Pm.dat")

Pp=0.0_dp   !Partikula karga berdnarekin
Pm=0.0_dp   !Partikula karga desberdinarekin
p=0.0_dp    !partikula kopurua

do n=1,1000


call equil(m)

do i=1,100
 k=i

 q=m(i,4)           !Zein da karga

 do j=1,99
     k=k+1
     if (k==101)then
          k=1
     end if


     d=abs(m(i,1)-m(k,1))    !Distantzia neurtu

       karg(j,1)=d                 !Sartu d
       karg(j,2)=m(k,4)            !Sartu bere karga

 end do

 do j=1,49   !tarteak
   d=real(j)*0.004_dp
   do k=1,99
     if (karg(k,1)>(d-0.004_dp).and.karg(k,1)<d)then      !Begiratu tartean badagoen
        if (karg(k,2)==q) then     !Sartu berdinen taldean
            Pp(j)=Pp(j)+1.0_dp
            p=p+1.0_dp
        else                   !Desberdinen taldean
           Pm(j)=Pm(j)+1.0_dp
           p=p+1.0_dp
        end if
     end if
   end do
 end do

end do

end do

s=0.0_dp
l=0.0_dp
do i=1,49
s=s+(Pp(i)/p)
l=l+(Pm(i)/p)
end do
print*,s
print*,l
print*,l+s

do i=1,49
  d=real(i)*0.004_dp
  s=Pp(i)+Pm(i)
  write(unit=11,fmt="(2f24.16)")d,Pp(i)/s
  write(unit=12,fmt="(2f24.16)")d,Pm(i)/s
end do



end program pant
