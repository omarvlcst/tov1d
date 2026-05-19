module ode_system
use constantes
implicit none
real(8) :: f, Q, r
contains

    subroutine f(t, y, dydt)
        use constantes
        implicit none
        real(8), intent(in) :: t
        real(8), intent(in) :: y(:), eps(:)
        real(8), intent(out) :: dydt(:)

        dydt(0)=4.d0*pi*eps(:)*t**2
        dydt(1)=-(eps(:)-y(1))*(y(0)+4*pi*t**3*y(1))/(t*(t-2*y(0)))
        dydt(2)=-y(2)*dydt(1)
    end subroutine f


    ! Crear antes las funciones a2 (g_rr), cs2(en funcion de las variables necesarias) y dadr (a es alpha)
    function Q(r,y,rho)
        use constantes
        implicit none
        real(8), intent(in) :: a2, cs2, dadr
        real(8) :: y(:), eps(:)
        real(8) :: r, Q1, Q2, Q3, Q

        Q1 = 4.d0 * pi * a2(r) * (5*eps + 9*y + (eps+y)/cs2)
        Q2 = -6* a2(r)/r**2
        Q3 = -4*dadr
        Q = Q1+Q2+Q3
    end function Q
        
end module ode_system