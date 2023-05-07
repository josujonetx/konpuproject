module funtzioak
use tipoak
      
public:: Plasma, Kanpoan        
            contains
            
        function Plasma(t, p) result(yprima)

        real(kind = dp), intent(in) :: t
        real(kind = dp), dimension(:), intent(in) :: p
        real(kind = dp), dimension(size(p)) :: yprima
        real(kind = dp), dimension(size(p)/6) :: x, y, z, vx, vy, vz, ax,ay, az,c
        integer:: i,j,n
        real(kind=dp):: r

        !posizioak eta abiadurak
        
        n=size(p)/6 !partikula kop.

        do i =1,n
                x(i)=p(6*i-5)
                y(i)=p(6*i-4)
                z(i)=p(6*i-3)
                vx(i)=p(6*i-2)
                vy(i)=p(6*i-1)
                vz(i)=p(6*i)
        enddo
        
        !kargen balioak
        
        do i=1,n/2
            c(i)=1.0_dp
            c(i+n/2)=-1.0_dp
        enddo
        
        do i=1,n
                do j=1, n
                     if (i/=j) then
                             r=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
                             ax(i)=ax(i)+c(i)*c(j)*(x(i)-x(j))/r**3
                             ay(i)=ay(i)+c(i)*c(j)*(y(i)-y(j))/r**3
                             az(i)=az(i)+c(i)*c(j)*(z(i)-z(j))/r**3
                     end if
                enddo
        enddo
        do i=1,n
                yprima(6*i-5)=vx(i)
                yprima(6*i-4)=vy(i)
                yprima(6*i-3)=vz(i)
                yprima(6*i-2)=ax(i)
                yprima(6*i-1)=ay(i)
                yprima(6*i)=az(i)
        enddo
        end function Plasma

          subroutine  Kanpoan(x,y,z,vx,vy,vz ) !Toroidetik atera ezin daitekenez
            real(kind=dp), intent(inout)::x,y,z, vx, vy, vz
            real(kind=dp):: r, theta, phy, x1,y1,z1
            r=sqrt((6-sqrt(x**2+y**2))**2+z**2)
            if (r>1.0_dp) then
                r=1.0_dp-r+1.0_dp !toroide barrura erreflejatu
                theta=atan(y/x)
                 if (x<0.0_dp) then
                         theta=theta+acos(-1.0)
                 end if
                phy=atan(z/(sqrt(x**2+y**2)-6))
                if (sqrt(x**2+y**2)<6.0) then
                        phy=phy+acos(-1.0)
                end if
                x=(6+r*cos(phy))*cos(theta)
                y=(6+r*cos(phy))*sin(theta)
                z=r*sin(phy)
                x1=-(6 - Sqrt(x**2 + y**2))*x/Sqrt(x**2 + y**2) !Toroidearen norma (x1,y1,z1)
                y1=-(6 - Sqrt(x**2 + y**2))*y/Sqrt(x**2 + y**2)
                z1=z
                x1=x1/sqrt(x1**2+y1**2+z1**2)
                y1=y1/sqrt(x1**2+y1**2+z1**2)
                z1=z1/(sqrt(x1**2+y1**2+z1**2))
                vx=vx-2*(vx*x1+vy*y1+vz*z1)*x1 !Proiekzioa eta talkaren ondoriozko aldaketa
                vy=vy-2*(vx*x1+vy*y1+vz*z1)*y1
                vz=vz-2*(vx*x1+vy*y1+vz*z1)*z1
            end if        

        end subroutine Kanpoan

end module funtzioak
