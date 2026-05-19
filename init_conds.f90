module boundary_conditions
    use constantes
    implicit none
    integer(4), parameter, public :: Num_pasos = 300000
    real(8), parameter, public :: r0  = 0.d0
    real(8), parameter, public :: r_max  = 30.d0
    real(8), parameter, public :: delta_r = (r_max - r0) / dble(Num_pasos-1)
    real(8), parameter, public :: p_thr  = 1.d-6
end module boundary_conditions