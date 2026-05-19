program main
    use constantes
    use grid
    use funciones
    use rk4
    implicit none

    integer :: i, n, Num_pasos, indice_estrella
    real(8) :: r, delta_r, r0, r_max, rho, epsilon, X_r, cs_sq, a0
    real(8) :: k2, lambda, alpha
    real(8) :: Radio_estrella, Masa_estrella, Presion_estrella, Lapso_estrella, X_estrella, constante_reescal
    real(8), allocatable :: y(:), y_out(:), y_estrella(:), r_out(:)

    n = 5                       ! numero de variables en el sistema de EDOs a integrar
    Num_pasos = 1000000
    r0 = 0.d0                   !punto inicial en coordenada r 
    r_max = 30.0d0
    delta_r = (r_max - r0) / dble(Num_pasos-1)
    call build_grid(Num_pasos,r0,delta_r,r_out)
    a0 = 0.5d0                  ! Parametro de regularidad cerca de r=0 para deformabilidad, hacer luego 0.1 y 1.0
    alpha= 1.d0
    allocate(y(n), y_out(n), y_estrella(n))

    ! Condiciones iniciales
    y(1) = alpha                      ! Alpha: función de lapso 
    y(2) = K_poly * rho_c0**(Gamma)   ! Presión central, Pin = K*rho**Gamma
    y(3) = 0.d0                       ! Masa inicial, para tomar el valor exacto en el origen r=0
    y(4) = a0*delta_r**2              ! Condicion inicial para H, H(0)
    y(5) = 2*a0*delta_r               ! Condicion inicial para beta, beta(0)

    write(*,*) 'Lapso central: ', y(1)
    write(*,*) 'Presion central: ', y(2)
    write(*,*) 'Masa central: ', y(3)
    write(*,*) 'Densidad de masa central: ', densidad_EOS(y(2))
    write(*,*) 'Velocidad del sonido central: ', dsqrt(cs2(y(2)))

    y_out=0.d0
    y_estrella=0.d0

    open(unit=10, file="tov.dat", status="replace")
    write(10,'(A)') "r[km]     m(r)     alpha(r)     P(r)       X(r)       rho(r)     cs^2       H(r)       Beta(r)"

    do i = 0, Num_pasos-1
        rho = densidad_EOS(y(2))
        X_r = func_X(r_out(i),y(3))                !dsqrt(r/(r-2.d0*y(3)))
        cs_sq = cs2(y(2))
        write(10,*) r_out(i)*const_geo, y(3), y(1), y(2), X_r, rho, cs_sq, y(4), y(5)
        call metodo_rk4(tov_rhs, r_out(i), y, delta_r, y_out, n)
        y = y_out   
        if (rho <= 1.d-6*rho_c0) exit
    end do

    indice_estrella = i
    Radio_estrella = r_out(i)                       !geometric units, multiply by M_km
    Lapso_estrella = y(1)
    Presion_estrella = y(2)
    Masa_estrella = y(3)
    X_estrella = dsqrt(Radio_estrella/(Radio_estrella-2*Masa_estrella))
    constante_reescal = dsqrt(1.d0-2*Masa_estrella/Radio_estrella)/Lapso_estrella

    write(*,*) 'La estrella se encontro en el paso', indice_estrella
    write(*,*) 'El lapso en el radio de la estrella es = ', Lapso_estrella
    write(*,*) 'La presion en el radio de la estrella es = ', Presion_estrella
    write(*,*) 'La densidad en el radio de la estrella es = ', densidad_EOS(Presion_estrella)
    write(*,*) 'La masa de la estrella es: M = ', Masa_estrella, 'Ms'
    write(*,*) 'El radio de la estrella en km es: R = ', Radio_estrella*const_geo
    write(*,*) 'La constante de reescalamiento es: ', constante_reescal

    k2 = tidal_love_number_k2(y(3),r_out(i),y(4),y(5),densidad_EOS(y(2)))
    lambda = tidal_deformability(k2,y(3),r_out(i))

    write(*,*) 'Tidal love number k2 = ', k2
    write(*,*) 'Tidal deformability Lambda =', lambda

    do i=indice_estrella,Num_pasos
        rho = densidad_EOS(y(2))
        X_r = func_X(r_out(i),y(3))
        cs_sq = cs2(y(2))
        write(10,*) r_out(i)*const_geo, y(3), 1.d0/(X_r*constante_reescal), y(2), X_r, rho, cs_sq, y(4), y(5)
        call metodo_rk4(tov_rhs, r_out(i), y, delta_r, y_out, n)
        y = y_out
    end do

    close(10)
    deallocate(r_out, y, y_out,y_estrella)
end program main