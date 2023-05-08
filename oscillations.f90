program oscillations

use tipoak
use egoerak
use m_spline_dp

real(kind=dp), dimension(1001):: bb
real(kind=dp), dimension(100,4) :: m
real(kind=dp),dimension(101) :: x,y,ys
real(kind=dp) :: V,t,ttot,Ez,Ep
integer :: i, j,k,l

bb(:)=0.0_dp
open(unit=11,status="replace",action="write",file="osc.dat")

do                                                               !100 sistema konparatu
   if (i>99) then
      exit
   end if

   call equil(m)                                                 !Oreka egoeratik abiatu
   m(:,2)=m(:,2)+m(:,4)*1.0_dp                                   !Ioiei bultzada bat eman

   ttot=0.0_dp

   x(1)=0.0_dp
   y(1)=2.0_dp
   do j=1,100                                                      !Sistema bakoitzean 101 puntu kalkulatu
      call ebol(m,0.04_dp,t)
      ttot=ttot+t
      x(j+1)=ttot
      y(j+1)=0.02_dp*sum(m(:,4)*m(:,2))                            !Abiaduren diferentzia
   end do

   Ez=0.5_dp*sum(m(:,2)**2)
   Ep=0.0_dp
   do k=1,100
   do l=1,100
      Ep=Ep+m(k,4)*m(l,4)*abs(m(k,1)-m(l,1))
   end do
   end do
   Ep=-0.5_dp*Ep

   if (Ez+Ep>3000) then
      !print*, Ez+Ep
      cycle
   end if

   !print*, "ondo"
   i=i+1

   call spline(x,y,101,ys)                                    !Batazbestekoa egiteko interpolazioa erabili
   do j=0,1000
      call splint(x,y,ys,101,0.004_dp*j,V)
      bb(j+1)=bb(j+1)+V
   end do
end do

bb(:)=bb(:)/100.0_dp                                         !Batazbestekoa egin

write(unit=11,fmt="(2f20.15)") 0.0_dp,2.0_dp

do j=0,1000
   write(unit=11, fmt="(2f20.15)") 0.004_dp*j, bb(j+1)
end do

end program oscillations
