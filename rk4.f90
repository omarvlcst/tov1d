module rk4
    implicit none
contains

    subroutine metodo_rk4(f, x, y, h, y_out, n)
        interface
            subroutine f(x, y, dydx, n)
                implicit none
                integer, intent(in) :: n
                real(8), intent(in) :: x
                real(8), intent(in) :: y(n)
                real(8), intent(out) :: dydx(n)
            end subroutine f
        end interface

        integer, intent(in) :: n
        real(8), intent(in) :: x, h
        real(8), intent(in) :: y(n)
        real(8), intent(out) :: y_out(n)

        real(8) :: k1(n), k2(n), k3(n), k4(n), y_temp(n)
        integer :: i

        call f(x, y, k1, n)

        do i = 1, n
            y_temp(i) = y(i) + 0.5d0 * h * k1(i)
        end do
        call f(x + 0.5d0 * h, y_temp, k2, n)

        do i = 1, n
            y_temp(i) = y(i) + 0.5d0 * h * k2(i)
        end do
        call f(x + 0.5d0 * h, y_temp, k3, n)

        do i = 1, n
            y_temp(i) = y(i) + h * k3(i)
        end do
        call f(x + h, y_temp, k4, n)

        do i = 1, n
            y_out(i) = y(i) + h * (k1(i) + 2.d0 * k2(i) + 2.d0 * k3(i) + k4(i)) / 6.d0
        end do
    end subroutine metodo_rk4

end module rk4