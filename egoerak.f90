module egoerak

use tipoak

public :: hasego, equil, sort_matrix, ebol 

contains

subroutine hasego(mo)

real(kind=dp), dimension(100) :: x,w,v,c
real(kind=dp), intent(out), dimension(:,:) :: mo
real(kind=dp) :: U0, V0, wp, vp, T
integer :: i,j

call random_seed()                !Zenbaki aleatorien sekuentzia aldatzeko ejekuzio bakoitzean

do i=1,50                         !Kargei balioak eman
   c(i)=1.0_dp
   c(i+50)=-1.0_dp
end do

U0=10.0_dp

do                                !Posizioak zehaztu
   call random_number(x)
   
   V0=0.0_dp
   do i=1,100
      do j=1,100
         V0=V0+c(i)*c(j)*abs(x(i)-x(j))
      end do
   end do
   V0=-0.5_dp*V0

   if (V0<=U0) then               !Energia potentziala energia osoa baino txikiagoa izan behar da
      exit
   end if
end do

vp=sqrt((U0-V0)/50)

call random_number(w)                 !Abiaduren esleipena
w=2*w-1

wp=0.0_dp
do i=1,100
   wp=wp+w(i)**2
end do
wp=sqrt(wp/100)

v=(vp/wp)*w

mo(:,1)=x                         !Partikula bakoitzaren informazioa bektore batean (x,v,c), egoera matrize batean (mo)
mo(:,2)=v
mo(:,3)=c

end subroutine hasego

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine equil(mo)                                  !Oreka egoera sortu 
    
real(kind=dp),dimension(100,3) :: m
real(kind=dp), dimension (100) :: a
real(kind=dp), dimension(:,:), intent(out) :: mo
real(kind=dp) :: ezk,esk, t
integer :: i,j

call hasego(m)                                   !Hasierako egoerak sortu

call sort_matrix(m)                              !Hasierako posizioak ordenatu

do i=1,100                                        !Azelerazioak kalkulatu
   ezk=0
   if (i>1) then
   do j=1,i-1 
      ezk=ezk+m(j,3)
   end do
   end if

   esk=0
   if (i<100) then
   do j=i+1,100
      esk=esk+m(j,3)
   end do
   end if

   a(i)=m(i,3)*(ezk-esk)
end do

mo(:,1:2)=m(:,1:2)                                    !mo (x,v,a,c) matrizea sortu
mo(:,3)=a(:)
mo(:,4)=m(:,3)

call ebol(mo,4.0_dp,t)

end subroutine equil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sort_matrix(mat)                                                          !matrizea ordenatzeko bubble algoritmoa
  real(kind=dp),intent(inout), dimension(:,:) :: mat
  integer :: i, j
  real(kind=dp), dimension(4) :: temp

  do i = 1, 100
    do j = 1, 99
      if (mat(j,1) > mat(j+1,1)) then
        temp = mat(j,:)
        mat(j,:) = mat(j+1,:)
        mat(j+1,:) = temp
      end if
    end do
  end do

end subroutine sort_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine ebol(mo,tmax,t)                                                    !denboraren eboluzioa

real(kind=dp), dimension(:,:), intent(inout) :: mo
real(kind=dp), intent(in) :: tmax
real(kind=dp), intent(out) :: t
real(kind=dp), dimension (4) :: aux
real(kind=dp) :: tau,plus,minus,a1,a2,barru
integer :: i,j,gert

t=0.0_dp
call sort_matrix(mo)

do  
   if (t>=tmax) then   
      exit
   end if

tau=10.0_dp
gert=-1

do i=1,99                                                                                !Partikulen arteko talken denborak
   if (mo(i+1,3)/=mo(i,3)) then
   barru=(mo(i+1,2)-mo(i,2))**2.0_dp-2.0_dp*(mo(i+1,1)-mo(i,1))*(mo(i+1,3)-mo(i,3))
   if (barru>=0.0_dp) then
      plus=(-(mo(i+1,2)-mo(i,2))+sqrt(barru))/(mo(i+1,3)-mo(i,3))
      if (plus<=tau .and. plus>1.0E-14_dp) then
         tau=plus
         gert=i 
      end if

      minus=(-(mo(i+1,2)-mo(i,2))-sqrt(barru))/(mo(i+1,3)-mo(i,3))
      if (minus<=tau .and. minus>1.0E-14_dp) then
         tau=minus
         gert=i
      end if
   end if
   else
   plus=-(mo(i+1,1)-mo(i,1))/(mo(i+1,2)-mo(i,2))
   if (plus<=tau .and. plus>1.0E-14_dp) then
      tau=plus
      gert=i
   end if
   end if
end do

barru=mo(1,2)**2.0_dp-2.0_dp*mo(1,1)*mo(1,3)                                                       !Ezker pareta
if (barru>=0.0_dp) then
   plus=(-mo(1,2)+sqrt(barru))/mo(1,3)
   if (plus<tau .and. plus>1.0E-14_dp) then
      tau=plus
      gert=0
   end if

   minus=(-mo(1,2)-sqrt(barru))/mo(1,3)
   if (minus<tau .and. minus>1.0E-14_dp) then
      tau=minus
      gert=0
   end if
end if

barru=mo(100,2)**2.0_dp-2.0_dp*(mo(100,1)-1)*mo(100,3)                                             !Eskuin pareta
if (barru>=0) then
   plus=(-mo(100,2)+sqrt(barru))/mo(100,3)
   if (plus<tau .and. plus>1.0E-14_dp) then
      tau=plus
      gert=100
   end if

   minus=(-mo(100,2)-sqrt(barru))/mo(100,3)
   if (minus<tau .and. minus>1.0E-14_dp) then
      tau=minus
      gert=100
   end if
end if

mo(:,1)=mo(:,1)+mo(:,2)*tau+mo(:,3)*(tau**2.0_dp)/2.0_dp                                            !Denboran eboluzioa
mo(:,2)=mo(:,2)+mo(:,3)*tau

if (gert>0 .and. gert<100) then                                                !Bi partikulen talka gertatzen bada
   a1=mo(gert,3)
   a2=mo(gert+1,3)
       
   aux=mo(gert,:) 
   mo(gert,:)=mo(gert+1,:)
   mo(gert+1,:)=aux
 
   mo(gert,3)=a2-2.0_dp*mo(gert,4)*mo(gert+1,4)
   mo(gert+1,3)=a1+2.0_dp*mo(gert,4)*mo(gert+1,4)
end if

if (gert==0) then                                                                      !Paretekin talka badago
   mo(1,2)=-mo(1,2)
end if

if (gert==100) then
   mo(100,2)=-mo(100,2)
end if

t=t+tau

end do

end subroutine ebol

end module egoerak
