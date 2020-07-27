#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      real_t  function getphirzfunc(radius)
      implicit none
      real_t radius
      getphirzfunc = radius*radius
      return
      end
      real_t  function getgradphirzfunc(radius)
      implicit none
      real_t radius
      getgradphirzfunc = two*radius
      return
      end
      real_t  function getlaplphirzfunc(radius)
      implicit none
      real_t radius
      getlaplphirzfunc = four
      return
      end
        subroutine GETPHI(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T probhi(0:0)
      integer iboxlo0
      integer iboxhi0
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getphipoint(phi(i),freq,x)
        
      enddo
        return
        end
        subroutine GETMAGRESIST(
     &           mag
     &           ,imaglo0
     &           ,imaghi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imaglo0
      integer imaghi0
      REAL_T mag(
     &           imaglo0:imaghi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      integer iboxlo0
      integer iboxhi0
      integer icomp
      integer whichmag
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getmagpointresist(mag(i),freq,x,
     $         icomp, whichmag)
        
      enddo
        return
        end
        subroutine GETMAGPOINTRESIST(
     &           mag
     &           ,freq
     &           ,xval
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T mag
      REAL_T freq(0:0)
      REAL_T xval(0:0)
      integer icomp
      integer whichmag
        REAL_T  x
        integer i
#if CH_SPACEDIM==2
        if(icomp.eq.2) then
           mag = zero
           return
        endif
#endif
        
        i = icomp
        if(whichmag.eq. 2) then
           
           x = freq(i)*xval(i)
#if CH_SPACEDIM==2
           mag = sin(y)
#elif CH_SPACEDIM==3
           mag = sin(y) + sin(z)
#else
           mag = x
#endif
        elseif(whichmag.eq. 1) then
           x = freq(icomp)*xval(icomp)
           mag = sin(x)
        elseif(whichmag.eq.0) then
           x = xval(icomp)
           mag = x*x
        elseif(whichmag.eq.4) then
           x = xval(icomp)
           if(icomp .eq. 0) then
              mag = x
           else
              mag = zero
           endif
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETETARESIST(
     &           eta
     &           ,ietalo0
     &           ,ietahi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,idir
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ietalo0
      integer ietahi0
      REAL_T eta(
     &           ietalo0:ietahi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      integer iboxlo0
      integer iboxhi0
      integer idir
      REAL_T eps
      integer whicheta
        integer i,jdir
        real_t x(0:CH_SPACEDIM-1)
        integer iv(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

           
           iv(0) = i
           do jdir = 0, CH_SPACEDIM-1
              if(idir .eq. jdir) then
                 x(jdir) = iv(jdir)*dx(jdir) + problo(jdir)
              else
                 x(jdir) = (iv(jdir)+half)*dx(jdir) + problo(jdir)
              endif
           enddo
          call getetapointresist(eta(i),freq,x,
     $         idir, eps, whicheta)
        
      enddo
        return
        end
        subroutine GETETAPOINTRESIST(
     &           eta
     &           ,freq
     &           ,xval
     &           ,idir
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T eta
      REAL_T freq(0:0)
      REAL_T xval(0:0)
      integer idir
      REAL_T eps
      integer whicheta
        REAL_T x
        if(whicheta.eq. 1) then
           
           x = freq(0)*xval(0)
           eta  = one + eps*(sin(x))
        elseif(whicheta.eq.0) then
           eta = one
        elseif(whicheta.eq.3) then
           eta = half
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETBETAVISCOUS(
     &           beta
     &           ,ibetalo0
     &           ,ibetahi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,eps
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,whichbeta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ibetalo0
      integer ibetahi0
      REAL_T beta(
     &           ibetalo0:ibetahi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T eps
      integer iboxlo0
      integer iboxhi0
      integer whichbeta
        integer i,jdir
        real_t x(0:CH_SPACEDIM-1)
        integer iv(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

           
           iv(0) = i
           do jdir = 0, CH_SPACEDIM-1
              x(jdir) = (iv(jdir)+half)*dx(jdir) + problo(jdir)
           enddo
          call getbetapointviscous(beta(i),freq,x, eps, whichbeta)
        
      enddo
        return
        end
        subroutine GETBETAPOINTVISCOUS(
     &           beta
     &           ,freq
     &           ,xval
     &           ,eps
     &           ,whichbeta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T beta
      REAL_T freq(0:0)
      REAL_T xval(0:0)
      REAL_T eps
      integer whichbeta
        REAL_T x
        if(whichbeta.eq. 1) then
           
           x = freq(0)*xval(0)
           beta  = one + eps*(sin(x))
        elseif(whichbeta.eq.0) then
           beta = one
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBRESIST(
     &           klb
     &           ,iklblo0
     &           ,iklbhi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0
      integer iklbhi0
      REAL_T klb(
     &           iklblo0:iklbhi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0
      integer iboxhi0
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getklbpointresist(klb(i),freq,x,
     $         alpha, beta, icomp, eps, whichmag, whicheta)
        
      enddo
        return
        end
        subroutine GETKLBPOINTRESIST(
     &           klb
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,beta
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T klb
      REAL_T freq(0:0)
      REAL_T xvec(0:0)
      REAL_T alpha
      REAL_T beta
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        REAL_T  fx
        REAL_T  x,  termone
        REAL_T  freqx, mag, divf, eta
        integer i
#if CH_SPACEDIM==2
        if(icomp.eq.2) then
           klb = zero
           return
        endif
#endif
        
        i = icomp
        call getetapointresist(eta,freq,xvec, icomp, eps, whicheta)
        call getmagpointresist(mag,freq,xvec, icomp, whichmag)
        
        freqx = freq(i)
        
        x = xvec(i)
        if((whichmag.eq. 2).and.(whicheta.eq.0)) then
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = -(freqy*freqy*sin(freqy*y))
#if CH_SPACEDIM==3
           divf = divf - (freqz*freqz*sin(freqz*z))
#endif
#endif
        elseif((whichmag.eq. 3).and.(whicheta.eq.0)) then
           call getlofphipoint(klb, freq, xvec, alpha, beta)
           goto  123
        elseif((whichmag.eq. 2).and.(whicheta.eq.1)) then
           
           fx = freqx*x
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = (freqy*cos(fy) - freqx*cos(fx))*(eps*freqy*cos(fy)) - freqy*freqy*eta*sin(fy)
#if CH_SPACEDIM==3
           divf = eps*freqy*cos(fy)*(freqy*cos(fy) - freqx*cos(fx))  - eta*freqy*freqy*sin(fy)
     $          + eps*freqz*cos(fz)*(freqz*cos(fz) - freqx*cos(fx))  - eta*freqz*freqz*sin(fz)
#endif
#endif
        elseif((whichmag.eq. 1).and.(whicheta.eq.1)) then
           termone =  
     $          freqx*cos(freqx*x)
          divf = eps*freqx*cos(freqx*x)*termone
     $         -freqx*freqx*sin(freqx*x)*eta
        elseif((whichmag.eq.0).and.(whicheta.eq.0)) then
           divf = two
        elseif((whichmag.eq.4).and.(whicheta.eq.0)) then
           divf = zero
        elseif((whichmag.eq.1).and.(whicheta.eq.0)) then
           divf = -freqx*freqx*sin(freqx*x)
        else
           call MayDay_Error()
        endif
        klb = alpha*mag  + beta*divf
  123   continue
        return
        end
        subroutine GETKLVVISCOUS(
     &           klb
     &           ,iklblo0
     &           ,iklbhi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,lambdafactor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0
      integer iklbhi0
      REAL_T klb(
     &           iklblo0:iklbhi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0
      integer iboxhi0
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T lambdafactor
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getklvpointviscous(klb(i),
     $         freq,x, alpha, beta, icomp, eps,
     $         whichvel, whicheta, whichlambda, lambdafactor)
        
      enddo
        return
        end
        subroutine GETKLVPOINTVISCOUS(
     &           klv
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,beta
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,lambdafactor
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T klv
      REAL_T freq(0:0)
      REAL_T xvec(0:0)
      REAL_T alpha
      REAL_T beta
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T lambdafactor
        REAL_T  x,    lambda
        REAL_T  freqx, vel, divf, eta
        real_t fx
        integer i
        
        i = icomp
        call getetapointresist(   eta,freq,xvec, icomp, eps, whicheta)
        if(whichlambda .eq. 2) then
           lambda = -lambdafactor*eta
        else
           call getetapointresist(lambda,freq,xvec, icomp, eps, whichlambda)
        endif
        call getmagpointresist(   vel,freq,xvec, icomp, whichvel)
        
        freqx = freq(i)
        
        x = xvec(i)
        
        fx = freqx*x
        if((whichvel.eq.2).and.(whicheta.eq.0)) then
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = -freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf -freqz*freqz*sin(fz)
#endif
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.3).and.(whichlambda.eq.3)) then
           divf = -three*half*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.0)) then
           divf = -three*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.2)) then
           divf = -(two - lambdafactor)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.1)) then
#if CH_SPACEDIM==1
           divf = one
#else
           divf =       -three*eta*freqx*freqx*sin(fx)
           divf = divf + three*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf +       eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf +       eps*freqx*freqz*cos(fx)*cos(fz)
#endif
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.2)) then
#if CH_SPACEDIM==1
           divf = one
#else
           divf =       -(two - lambdafactor)*eta*freqx*freqx*sin(fx)
           divf = divf + (two - lambdafactor)*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf -       (lambdafactor)*eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf -       (lambdafactor)*eps*freqx*freqz*cos(fx)*cos(fz)
#endif
#endif
        else if((whichvel.eq.2).and.(whicheta.eq.1)) then
#if  CH_SPACEDIM==1
           divf  = one
#else
           divf =       - eta*freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf  - eta*freqz*freqz*sin(fz)
#endif
#endif
#if  CH_SPACEDIM==1
           divf = one
#else
           divf = divf + eps*freqy*cos(fy)*(freqy*cos(fy) + freqx*cos(fx))
#if CH_SPACEDIM==3
           divf = divf + eps*freqz*cos(fz)*(freqz*cos(fz) + freqx*cos(fx))
#endif
#endif
        else if(whichvel.eq.4) then
           divf = zero
        else
           call MayDay_Error()
        endif
        klv = alpha*vel  + beta*divf
        return
        end
        subroutine GETPHIPOINT(
     &           phi
     &           ,freq
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T freq(0:0)
      REAL_T x(0:0)
        phi = sin(freq(0)*x(0))
        return
        end
        subroutine GETLOFPHIZPOLY(
     &           lofphi
     &           ,x
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T x(0:0)
      REAL_T alpha
      REAL_T beta
        real_t phi, laplphi
        real_t dist
        external getlaplphirzfunc
        real_t getlaplphirzfunc
        external getphirzfunc
        real_t getphirzfunc
        dist = abs(x(0))
        phi = getphirzfunc(dist)
        laplphi = getlaplphirzfunc(dist)
        lofphi = alpha*phi + beta*laplphi
        return
        end
        subroutine GETPHIRZPOLY(
     &           phi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T phi
      REAL_T x(0:0)
        real_t dist
        external getphirzfunc
        real_t getphirzfunc
        dist =abs(x(0))
        phi = getphirzfunc(dist)
        return
        end
      subroutine GETGRADPHIRZPOLY(
     &           gradphi
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:0)
      REAL_T x(0:0)
        real_t dist
        external getgradphirzfunc
        real_t getgradphirzfunc
        dist = abs(x(0))
        
        gradphi(0) = getgradphirzfunc(dist)
        return
        end
        subroutine GETGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:0)
      REAL_T freq(0:0)
      REAL_T x(0:0)
        
        gradphi(0) = freq(0) * cos(freq(0)*x(0))                                       
        return
        end
        subroutine GETLOFPHI(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,beta
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0
      integer ilofphihi0
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T probhi(0:0)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0
      integer iboxhi0
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getlofphipoint(lofphi(i),freq,x,alpha,beta)
        
      enddo
        return
        end
        subroutine GETLOFPHIPOINT(
     &           lofphi
     &           ,freq
     &           ,x
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T lofphi
      REAL_T freq(0:0)
      REAL_T x(0:0)
      REAL_T alpha
      REAL_T beta
        real_t fac,phi
        fac = -(freq(0)**2)
        phi = (sin(freq(0)*x(0)))
        lofphi = fac*phi
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETDBGPHI(
     &           dbgphi
     &           ,idbgphilo0
     &           ,idbgphihi0
     &           ,beta
     &           ,ibetalo0
     &           ,ibetahi0
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idbgphilo0
      integer idbgphihi0
      REAL_T dbgphi(
     &           idbgphilo0:idbgphihi0)
      integer ibetalo0
      integer ibetahi0
      REAL_T beta(
     &           ibetalo0:ibetahi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T probhi(0:0)
      REAL_T alpha
      integer iboxlo0
      integer iboxhi0
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getdbgphipoint(dbgphi(i),
     &        beta(i),freq,x,alpha)
        
      enddo
        return
        end
        subroutine GETDBGPHIPOINT(
     &           dbgphi
     &           ,beta
     &           ,freq
     &           ,x
     &           ,alpha
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dbgphi
      REAL_T beta
      REAL_T freq(0:0)
      REAL_T x(0:0)
      REAL_T alpha
        real_t gradphi(0:CH_SPACEDIM-1),gradbeta(0:CH_SPACEDIM-1)
        real_t alphaphiplusbetalapphi,gradbetadotgradphi
        call getbetapoint(beta,freq,x)
        call getlofphipoint(alphaphiplusbetalapphi,freq,x,alpha,beta)
        call getgradbetapoint(gradbeta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        gradbetadotgradphi = gradbeta(0)*gradphi(0)
        dbgphi = alphaphiplusbetalapphi
        dbgphi = dbgphi + gradbetadotgradphi
        return
        end
        subroutine GETBETAPOINT(
     &           beta
     &           ,freq
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T beta
      REAL_T freq(0:0)
      REAL_T x(0:0)
        beta = x(0)*x(0)
        return
        end
        subroutine GETGRADBETAPOINT(
     &           gradbeta
     &           ,freq
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradbeta(0:0)
      REAL_T freq(0:0)
      REAL_T x(0:0)
        integer idir
        do idir = 0, CH_SPACEDIM-1
            gradbeta(idir) = two*x(idir)
        enddo
        return
        end
        subroutine GETBETAGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T gradphi(0:0)
      REAL_T freq(0:0)
      REAL_T x(0:0)
        integer idir
        real_t beta
        call getbetapoint(beta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETSRC(
     &           src
     &           ,isrclo0
     &           ,isrchi0
     &           ,freq
     &           ,dx
     &           ,diffconst
     &           ,problo
     &           ,probhi
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0
      integer isrchi0
      REAL_T src(
     &           isrclo0:isrchi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T diffconst
      REAL_T problo(0:0)
      REAL_T probhi(0:0)
      integer iboxlo0
      integer iboxhi0
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          call getsrcpoint(src(i),freq,x,diffconst)
        
      enddo
        return
        end
        subroutine GETSRCPOINT(
     &           src
     &           ,freq
     &           ,x
     &           ,diffconst
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T src
      REAL_T freq(0:0)
      REAL_T x(0:0)
      REAL_T diffconst
        real_t fac,phi
        fac = -(freq(0)**2)
        phi = (sin(freq(0)*x(0)))
        src = (-fac*diffconst)*phi
        return
        end
