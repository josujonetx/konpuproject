module rk4
public ::  rk4_paso_dp

contains

subroutine rk4_paso_dp(x,y,f,h)
use tipoak

real(kind=dp), intent(inout)                :: x
real(kind=dp), dimension(:), intent(inout)  :: y
real(kind=dp), intent(in)                   :: h
interface
   function f(x,y) result(yprima)
    use tipoak
    real(kind=dp), intent(in)               :: x
    real(kind=dp), dimension(:), intent(in) :: y
    real(kind=dp), dimension(size(y))       :: yprima
   end function f
end interface
real(kind=dp), dimension(size(y))           :: k1, k2, k3, k4

   k1 = h * f(x,y)                        ! 
   k2 = h * f(x+0.5_dp*h , y+0.5_dp*k1)   ! 
   k3 = h * f(x+0.5_dp*h , y+0.5_dp*k2)   ! 
   k4 = h * f(x+h , y+k3)                 ! 

   x = x + h
   y = y + (k1 + 2*k2 + 2*k3 +k4) / 6.0_dp    ! Salto definitivo

end subroutine rk4_paso_dp

end module rk4
