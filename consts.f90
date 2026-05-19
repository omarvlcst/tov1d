module constantes
    implicit none
    real(8), parameter, public :: pi        = 4.d0*datan(1.d0)
    real(8), parameter, public :: G         = 6.6743015d-20    ! constante de gravitacion universal en km³/(kg*s²)
    real(8), parameter, public :: c         = 299792.458d0     ! velocidad de la luz en km/s
    real(8), parameter, public :: M_sun     = 1.98847d30       ! masa del sol en kg 
    real(8), parameter, public :: M_km      = G*M_sun/c**2     ! una masa solar convertida a unidades geometricas
    real(8), parameter, public :: R_km      = 14.16d0          ! radio de la NS en km
    real(8), parameter, public :: M_star    = 1.4d0 * M_km     ! masa de la NS en km (unidades de codigo)
    real(8), parameter, public :: K_poly    = 100.d0           ! constante politropica
    real(8), parameter, public :: Gamma     = 2.d0             ! indice politropico        
    real(8), parameter, public :: rho_c0    = 1.28d-3          ! densidad central inicial de la NS en unidades de codigo 
    real(8), parameter, public :: const_geo = M_km             ! constante para reescalar las distancias de unidades geom a cgs (Rezzolla, 2012)
end module constantes