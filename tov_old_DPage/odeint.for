      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad)
       INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX  
       REAL*8 eps,h1,hmin,x1,x2,TINY    
       PARAMETER (MAXSTP=100000,NMAX=50,KMAXX=100000,TINY=1.d-30)  
       INTEGER i,kmax,kount,nstp  
       REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),  
     X yp(NMAX,KMAXX),yscal(NMAX),ystart(nvar)
       COMMON /path/ kmax,kount,dxsav,xp,yp  
        x=x1  
        h=sign(h1,x2-x1)  
        nok=0  
        nbad=0  
        kount=0  
        do i=1,nvar  
         y(i)=ystart(i) 
        enddo 
        if (kmax.gt.0) xsav=x-2.d0*dxsav  
        do nstp=1,MAXSTP  
         call derivs(x,y,dydx)
         do i=1,nvar  
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY  
	       enddo   
         if(kmax.gt.0)then  
           if(abs(x-xsav).gt.abs(dxsav)) then  
             if(kount.lt.kmax-1)then  
               kount=kount+1  
               xp(kount)=x  
               do i=1,nvar  
                 yp(i,kount)=y(i)  
	             enddo   
               xsav=x 
             endif  
           endif  
         endif  
         if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x  
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext)  
         if(hdid.eq.h)then  
           nok=nok+1  
         else  
           nbad=nbad+1  
         endif  
         if((x-x2)*(x2-x1).ge.0.d0)then  
           do i=1,nvar  
             ystart(i)=y(i)  
	         enddo  
           if(kmax.ne.0)then  
             kount=kount+1  
             xp(kount)=x 
             do i=1,nvar  
               yp(i,kount)=y(i)
	           enddo  
           endif  
           return  
         endif  
         if(abs(hnext).lt.hmin) pause  
     X	     'stepsize smaller than minimum in odeint'  
         h=hnext  
        enddo  
       pause 'too many steps in odeint'  
       return  
       END