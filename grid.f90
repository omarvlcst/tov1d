module grid
    use constantes
    implicit none
contains

    subroutine build_grid(N_r, r_ini, dr, r)
        implicit none
        integer(4)             :: i
        integer(4), intent(in) :: N_r
        real(8), intent(in)               :: dr, r_ini
        real(8), intent(out), allocatable :: r(:)

        allocate(r(0:N_r))

        do i = 0, N_r
            r(i) = r_ini + i * dr
        end do
    end subroutine build_grid

end module grid