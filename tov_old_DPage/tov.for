      PROGRAM TOV
cf2py intent(in) w0
      integer nvar,nok,nbad,NMAX,KMAXX,kmax,kount,k,nrhs
      PARAMETER (NMAX=50,KMAXX=100000,nvar=2)
      real*8 dxsav, xp(KMAXX), yp(NMAX,KMAXX)
      real*8 x1,x2,eps,h1,hmin,pi,ystart(nvar),w0
      PARAMETER (w0=3.034d0)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      COMMON nrhs
      EXTERNAL derivs, rkqs, odeint
      nrhs=0
      dxsav = 1.d-12
      kmax = KMAXX
      ystart(1) = w0
      ystart(2) = 0.0d0
      x1 = 0.0001d0
      x2 = 1.d0
      eps = 1.d-13
      h1=1.d-10
      hmin=0.d0
      call odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad)
      do k=1,kount
        if(yp(1,k).gt.1.d0) write(*,'(1x,i4,1x,f10.4,2x,2f14.6)') 
     X                       k, xp(k), yp(1,k), xp(k)/yp(2,k)
      end do
      END PROGRAM TOV

      SUBROUTINE derivs(x,y,dydx)
       integer NMAX, nrhs
       real pi,fw1,fw2
       PARAMETER (NMAX=50, pi=4.d0*datan(1.d0))
       real*8 x,y(NMAX),dydx(NMAX)
       COMMON nrhs
       nrhs=nrhs+1
       fw1 = y(2)+4.d0*pi*x**3*(y(1)-1.d0)
       fw2 = x*(x-2.d0*y(2))
       dydx(1) = -(2.d0*y(1)-1.d0)*fw1/fw2
       dydx(2) = 4.d0*pi*x**2*y(1)
       RETURN
       END SUBROUTINE

       INCLUDE 'odeint.for'
       INCLUDE 'rkck.for'
       INCLUDE 'rkqs.for' 