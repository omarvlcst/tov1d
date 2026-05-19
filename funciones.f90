module funciones
use constantes
implicit none
contains

    subroutine tov_rhs(r, y, dydr, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: r
        real(8), intent(in) :: y(n)
        real(8), intent(out) :: dydr(n)
        real(8) :: alpha, P, m, rho, epsilon, H_r, beta_r, dalphadr, dPdr, dmdr, dHdr, dbetadr
        real(8) :: denom, Qr, Sr, cs2, g_rr

        alpha = y(1)        ! g_tt
        P = y(2)            ! P=P(r) presión
        m = y(3)            ! m=m(r) masa encerrada
        H_r = y(4)          ! H=H(r)
        beta_r = y(5)       ! beta=beta(r)

        if (P <= 0.d0 .or. r <= 0.d0) then
            dydr(1) = 0.d0
            dydr(2) = 0.d0
            dydr(3) = 0.d0
            dydr(4) = 0.d0
            dydr(5) = 0.d0
            return
        end if

        rho = (P / K_poly)**(1.d0 / Gamma)      ! Densidad de masa EOS
        epsilon = rho + K_poly* rho**(Gamma)/(Gamma-1.d0)          ! Densidad de energía (e)
        cs2 = K_poly * Gamma * rho**(Gamma-1.d0)         ! Velocidad del sonido relativista
        denom = r * (r - 2.d0 * m)              ! Denominador en Ecuacion para P y alpha
        g_rr = 1.d0/(1.d0-2.d0*m/r)                    ! g_rr = 1 - 2m/r
        Qr = -2.d0 * pi * (5.d0 * epsilon + 9.d0 * P + (epsilon + P) / cs2)  &
        + 3.d0 / r**2                                                     &
        + 2.d0 * g_rr * ((m / r**2 + 4.d0 * pi * r * P)**2)
        Sr = 2.d0 * beta_r / r * g_rr * (-1.d0 + m/r + 2.d0 * pi * r**2 * (epsilon-P))
    
        if (denom <= 1.d-3) then
            dydr(1) = 0.d0
            dydr(2) = 0.d0
            dydr(3) = 0.d0
            dydr(4) = 0.d0
            dydr(5) = 0.d0
            return
        end if

        if (r <= 0.03d0) then
            dalphadr = alpha * (4.d0 * pi * epsilon * r + 12.d0 * pi * r * P) / (3.d0- 8.d0*pi * epsilon * r**2)
            dPdr = -(epsilon + P) * (4.d0 * pi * epsilon * r + 12.d0 * pi * r * P) / (3.d0- 8.d0*pi * epsilon * r**2)
        else 
            dalphadr = alpha * (m + 4.d0 * pi * r**3 * P) / denom
            dPdr = -(epsilon + P) * (m + 4.d0 * pi * r**3 * P) / denom
        end if

        dmdr = 4.d0 * pi * r**2 * epsilon
        dHdr = beta_r
        dbetadr = 2.d0 * g_rr * H_r * Qr + Sr

        dydr(1) = dalphadr
        dydr(2) = dPdr
        dydr(3) = dmdr
        dydr(4) = dHdr
        dydr(5) = dbetadr

    end subroutine tov_rhs

    real(8) function func_X(radio,masa) result(X)
        implicit none
        real(8), intent(in) :: radio, masa
        if (radio==0.d0.and.masa==0.d0) then
            X = 1.d0
        else
            X = dsqrt(radio/(radio-2.d0*masa))
        end if
    end function func_X
    
    real(8) function densidad_EOS(P) result(rho)
        implicit none
        real(8), intent(in) :: P
        if (P > 0.d0) then
            rho = (P / K_poly)**(1.d0 / Gamma)
        else
            rho = 0.d0
        end if
    end function densidad_EOS

    real(8) function densidad_de_energia(P) result(epsilon)
        implicit none
        real(8), intent(in) :: P
        if (P > 0.d0) then
            epsilon = densidad_EOS(P) + P/(Gamma-1.d0)
        else 
            epsilon = 0.d0
        end if
    end function densidad_de_energia
    
    real(8) function cs2(P) result(cs_2)
        implicit none
        real(8), intent(in) :: P
        if (P > 0.d0) then
            !cs_2 = (Gamma * P) / (P+densidad_de_energia(P))
            cs_2 = K_poly * Gamma * densidad_EOS(P)**(Gamma-1.d0)
        else 
            cs_2 = 0.d0
        end if
    end function cs2

    real(8) function tidal_love_number_k2(m,r,H,beta,rho) result(k2)
        implicit none
        real(8), intent(in) :: m,r,H,beta,rho
        real(8) :: compactness, y, k2_1, k2_2, k2_3, k2_4, k2_5
        compactness = m/r 
        y = r * beta / H
        k2_1 = (8.d0 * compactness**5)/5 *(1.d0-2.d0*compactness)**2
        k2_2 = 2.d0+2.d0*compactness*(y-1.d0)-y
        k2_3 = 2.d0 * compactness*(6.d0-3.d0*y+3.d0*compactness*(5.d0*y-8.d0))
        k2_4 = 4.d0 * compactness**3 * (13.d0 - 11.d0*y+ compactness * (3.d0*y-2.d0) &
        +2.d0*compactness**2 *(1.d0+y))
        k2_5 = 3.d0 * (1.d0-2.d0*compactness)**2 * (2.d0 - y + 2.d0 * compactness * (y-1.d0))* dlog(1.d0-2.d0*compactness)
        k2 = k2_1*k2_2/(k2_3+k2_4+k2_5)
    end function tidal_love_number_k2

    real(8) function tidal_deformability(lnum_k2,m,r) result(lambda)
        implicit none
        real(8), intent(in) :: lnum_k2, m, r
        lambda=2.d0/3.d0*lnum_k2*(r/m)**5        
    end function tidal_deformability

end module funciones