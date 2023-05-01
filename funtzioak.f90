module funtzioak
use tipoak
      public:: Plasma, Potentziala, Deribatux, Deribatuy, Deribatuz 
        
            contains
            
        function Plasma(t, p) result(yprima)

        real(kind = dp), intent(in) :: t
        real(kind = dp), dimension(:), intent(in) :: p
        real(kind = dp), dimension(size(p)) :: yprima
        real(kind = dp), dimension(size(p)/6) :: x, y, z, vx, vy, vz, ax,ay, az,c
        integer:: i,j
        
        !posizioak eta abiadurak
        
        do i =1,size(y)/6
                x(i)=p(6*i-5)
                y(i)=p(6*i-4)
                z(i)=p(6*i-3)
                vx(i)=p(6*i-2)
                vy(i)=p(6*i-1)
                vz(i)=p(6*i)
        enddo
        
        !kargen balioak
        
        do i=1,size(p)/12
            c(i)=1
            c(i+25)=-1
        enddo
        
        do i=1,size(p)/6
                do j=1, size(p)/6
                     if (i/=j) then
                             ax(i)=ax(i)+c(i)*c(j)*(x(i)-x(j))/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2))**3
                             ay(i)=ay(i)+c(i)*c(j)*(y(i)-y(j))/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2))**3
                             az(i)=az(i)+c(i)*c(j)*(z(i)-z(j))/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2))**3
                     end if
                enddo
        enddo
        do i=1,size(y)/6
                yprima(6*i-5)=vx(i)
                yprima(6*i-4)=vy(i)
                yprima(6*i-3)=vz(i)
                yprima(6*i-2)=ax(i)-Deribatux(x(i),y(i),z(i), Potentziala)
                yprima(6*i-1)=ay(i)-Deribatuy(x(i),y(i),z(i), Potentziala)
                yprima(6*i)=az(i)-Deribatuz(x(i),y(i),z(i), Potentziala)
        enddo
        end function Plasma
        
        function Potentziala(x,y,z)
            real(kind=dp), intent(in)::x,y,z
            real(kind=dp):: Potentziala
        
            Potentziala=Sqrt[(6 - Sqrt[x^2 + y^2])^2 + z^2] !Donut bat
        
        end function Potentziala
        
        function Deribatux(x,y,z,f)
        
        real(kind=dp), intent(in)::x,y,z
        real(kind=dp):: Deribatux
        interface
            function f(x,y,z)
                real(kind=dp), intent(in):: x,y,z
                real(kind=dp):: f
            end function f
        end interface
        h=1.0E-7
        Deribatux=(f(x+h,y,z)-f(x-h,y,z))/2/h
        
        end function Deribatux
        
        function Deribatuy(x,y,z,f)
        
        real(kind=dp), intent(in)::x,y,z
        real(kind=dp):: Deribatuy
        interface
            function f(x,y,z)
                real(kind=dp), intent(in):: x,y,z
                real(kind=dp):: f
            end function f
        end interface
        
        h=1.0E-7
        Deribatuy=(f(x,y+h,z)-f(x,y-h,z))/2/h
        
        end function Deribatuy
        
        function Deribatuz(x,y,z,f)
        
        real(kind=dp), intent(in)::x,y,z
        real(kind=dp):: Deribatux
        interface
            function f(x,y,z)
                real(kind=dp), intent(in):: x,y,z
                real(kind=dp):: f
            end function f
        end interface
        
        h=1.0E-7
        Deribatux=(f(x,y,z+h)-f(x,y,z+h))/2/h
        
        end function Deribatuz
end module funtzioak
