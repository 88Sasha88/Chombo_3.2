      REAL*8  function getphirzfunc(radius)
      implicit none
      REAL*8 radius
      getphirzfunc = radius*radius
      return
      end
      REAL*8  function getgradphirzfunc(radius)
      implicit none
      REAL*8 radius
      getgradphirzfunc = (2.0d0)*radius
      return
      end
      REAL*8  function getlaplphirzfunc(radius)
      implicit none
      REAL*8 radius
      getlaplphirzfunc = (4.0d0)
      return
      end
        subroutine GETPHI(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getphipoint(phi(i,j),freq,x)
      enddo
      enddo
        return
        end
        subroutine GETMAGRESIST(
     &           mag
     &           ,imaglo0,imaglo1
     &           ,imaghi0,imaghi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,icomp
     &           ,whichmag
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer imaglo0,imaglo1
      integer imaghi0,imaghi1
      REAL*8 mag(
     &           imaglo0:imaghi0,
     &           imaglo1:imaghi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      integer whichmag
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getmagpointresist(mag(i,j),freq,x,
     $         icomp, whichmag)
      enddo
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
      REAL*8 mag
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      integer icomp
      integer whichmag
        REAL*8  x,y
        integer i,j
        if(icomp.eq.2) then
           mag = (0.0d0)
           return
        endif
        i = icomp
        j = max(1-icomp, 0)
        if(whichmag.eq. 2) then
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
           mag = sin(y)
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
              mag = (0.0d0)
           endif
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETETARESIST(
     &           eta
     &           ,ietalo0,ietalo1
     &           ,ietahi0,ietahi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,idir
     &           ,eps
     &           ,whicheta
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ietalo0,ietalo1
      integer ietahi0,ietahi1
      REAL*8 eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idir
      REAL*8 eps
      integer whicheta
        integer i,j,jdir
        REAL*8 x(0:2 -1)
        integer iv(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
           iv(0) = i
           iv(1) = j
           do jdir = 0, 2 -1
              if(idir .eq. jdir) then
                 x(jdir) = iv(jdir)*dx(jdir) + problo(jdir)
              else
                 x(jdir) = (iv(jdir)+(0.500d0))*dx(jdir) + problo(jdir)
              endif
           enddo
          call getetapointresist(eta(i,j),freq,x,
     $         idir, eps, whicheta)
      enddo
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
      REAL*8 eta
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      integer idir
      REAL*8 eps
      integer whicheta
        REAL*8 x, y
        if(whicheta.eq. 1) then
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           eta  = (1.0d0) + eps*(sin(x) + sin(y))
        elseif(whicheta.eq.0) then
           eta = (1.0d0)
        elseif(whicheta.eq.3) then
           eta = (0.500d0)
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETBETAVISCOUS(
     &           beta
     &           ,ibetalo0,ibetalo1
     &           ,ibetahi0,ibetahi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,eps
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,whichbeta
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL*8 beta(
     &           ibetalo0:ibetahi0,
     &           ibetalo1:ibetahi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 eps
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer whichbeta
        integer i,j,jdir
        REAL*8 x(0:2 -1)
        integer iv(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
           iv(0) = i
           iv(1) = j
           do jdir = 0, 2 -1
              x(jdir) = (iv(jdir)+(0.500d0))*dx(jdir) + problo(jdir)
           enddo
          call getbetapointviscous(beta(i,j),freq,x, eps, whichbeta)
      enddo
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
      REAL*8 beta
      REAL*8 freq(0:1)
      REAL*8 xval(0:1)
      REAL*8 eps
      integer whichbeta
        REAL*8 x, y
        if(whichbeta.eq. 1) then
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           beta  = (1.0d0) + eps*(sin(x) + sin(y))
        elseif(whichbeta.eq.0) then
           beta = (1.0d0)
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBRESIST(
     &           klb
     &           ,iklblo0,iklblo1
     &           ,iklbhi0,iklbhi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL*8 klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL*8 eps
      integer whichmag
      integer whicheta
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getklbpointresist(klb(i,j),freq,x,
     $         alpha, beta, icomp, eps, whichmag, whicheta)
      enddo
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
      REAL*8 klb
      REAL*8 freq(0:1)
      REAL*8 xvec(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer icomp
      REAL*8 eps
      integer whichmag
      integer whicheta
        REAL*8  fx,fy
        REAL*8  x,y,  termone
        REAL*8  freqx,freqy, mag, divf, eta
        integer i,j
        if(icomp.eq.2) then
           klb = (0.0d0)
           return
        endif
        i = icomp
        j = max(1-icomp, 0)
        call getetapointresist(eta,freq,xvec, icomp, eps, whicheta)
        call getmagpointresist(mag,freq,xvec, icomp, whichmag)
        freqx = freq(i)
        freqy = freq(j)
        x = xvec(i)
        y = xvec(j)
        if((whichmag.eq. 2).and.(whicheta.eq.0)) then
           divf = -(freqy*freqy*sin(freqy*y))
        elseif((whichmag.eq. 3).and.(whicheta.eq.0)) then
           call getlofphipoint(klb, freq, xvec, alpha, beta)
           goto  123
        elseif((whichmag.eq. 2).and.(whicheta.eq.1)) then
           fx = freqx*x
           fy = freqy*y
           divf = (freqy*cos(fy) - freqx*cos(fx))*(eps*freqy*cos(fy)) - 
     &freqy*freqy*eta*sin(fy)
        elseif((whichmag.eq. 1).and.(whicheta.eq.1)) then
           termone =  
     $          freqx*cos(freqx*x) +
     $          freqy*cos(freqy*y)
          divf = eps*freqx*cos(freqx*x)*termone
     $         -freqx*freqx*sin(freqx*x)*eta
        elseif((whichmag.eq.0).and.(whicheta.eq.0)) then
           divf = (2.0d0)
        elseif((whichmag.eq.4).and.(whicheta.eq.0)) then
           divf = (0.0d0)
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
     &           ,iklblo0,iklblo1
     &           ,iklbhi0,iklbhi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
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
      integer iklblo0,iklblo1
      integer iklbhi0,iklbhi1
      REAL*8 klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer icomp
      REAL*8 eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL*8 lambdafactor
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getklvpointviscous(klb(i,j),
     $         freq,x, alpha, beta, icomp, eps,
     $         whichvel, whicheta, whichlambda, lambdafactor)
      enddo
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
      REAL*8 klv
      REAL*8 freq(0:1)
      REAL*8 xvec(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer icomp
      REAL*8 eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL*8 lambdafactor
        REAL*8  x,y,    lambda
        REAL*8  freqx,freqy, vel, divf, eta
        REAL*8 fx,fy
        integer i,j
        i = icomp
        j = max(1-icomp, 0)
        call getetapointresist(   eta,freq,xvec, icomp, eps, whicheta)
        if(whichlambda .eq. 2) then
           lambda = -lambdafactor*eta
        else
           call getetapointresist(lambda,freq,xvec, icomp, eps, whichlam
     &bda)
        endif
        call getmagpointresist(   vel,freq,xvec, icomp, whichvel)
        freqx = freq(i)
        freqy = freq(j)
        x = xvec(i)
        y = xvec(j)
        fx = freqx*x
        fy = freqy*y
        if((whichvel.eq.2).and.(whicheta.eq.0)) then
           divf = -freqy*freqy*sin(fy)
        else if((whichvel.eq.1).and.(whicheta.eq.3).and.(whichlambda.eq.
     &3)) then
           divf = -(3.0d0)*(0.500d0)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.
     &0)) then
           divf = -(3.0d0)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.
     &2)) then
           divf = -((2.0d0) - lambdafactor)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.
     &1)) then
           divf =       -(3.0d0)*eta*freqx*freqx*sin(fx)
           divf = divf + (3.0d0)*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf +       eps*freqx*freqy*cos(fx)*cos(fy)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.
     &2)) then
           divf =       -((2.0d0) - lambdafactor)*eta*freqx*freqx*sin(fx
     &)
           divf = divf + ((2.0d0) - lambdafactor)*eps*freqx*freqx*cos(fx
     &)*cos(fx)
           divf = divf -       (lambdafactor)*eps*freqx*freqy*cos(fx)*co
     &s(fy)
        else if((whichvel.eq.2).and.(whicheta.eq.1)) then
           divf =       - eta*freqy*freqy*sin(fy)
           divf = divf + eps*freqy*cos(fy)*(freqy*cos(fy) + freqx*cos(fx
     &))
        else if(whichvel.eq.4) then
           divf = (0.0d0)
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
      REAL*8 phi
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
        phi = sin(freq(0)*x(0))
     &                * sin(freq(1)*x(1))
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
      REAL*8 lofphi
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 beta
        REAL*8 phi, laplphi
        REAL*8 dist
        external getlaplphirzfunc
        REAL*8 getlaplphirzfunc
        external getphirzfunc
        REAL*8 getphirzfunc
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
      REAL*8 phi
      REAL*8 x(0:1)
        REAL*8 dist
        external getphirzfunc
        REAL*8 getphirzfunc
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
      REAL*8 gradphi(0:1)
      REAL*8 x(0:1)
        REAL*8 dist
        external getgradphirzfunc
        REAL*8 getgradphirzfunc
        dist = abs(x(0))
        gradphi(0) = getgradphirzfunc(dist)
        gradphi(1) = (0.0d0)
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
      REAL*8 gradphi(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
        gradphi(0) = freq(0) * cos(freq(0)*x(0)) * sin(freq(1)*x(1))
        gradphi(1) = freq(1) * sin(freq(0)*x(0)) * cos(freq(1)*x(1))    
        return
        end
        subroutine GETLOFPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0,ilofphilo1
      integer ilofphihi0,ilofphihi1
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      REAL*8 alpha
      REAL*8 beta
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getlofphipoint(lofphi(i,j),freq,x,alpha,beta)
      enddo
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
      REAL*8 lofphi
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 alpha
      REAL*8 beta
        REAL*8 fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2)
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1)))
        lofphi = fac*phi
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETDBGPHI(
     &           dbgphi
     &           ,idbgphilo0,idbgphilo1
     &           ,idbgphihi0,idbgphihi1
     &           ,beta
     &           ,ibetalo0,ibetalo1
     &           ,ibetahi0,ibetahi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer idbgphilo0,idbgphilo1
      integer idbgphihi0,idbgphihi1
      REAL*8 dbgphi(
     &           idbgphilo0:idbgphihi0,
     &           idbgphilo1:idbgphihi1)
      integer ibetalo0,ibetalo1
      integer ibetahi0,ibetahi1
      REAL*8 beta(
     &           ibetalo0:ibetahi0,
     &           ibetalo1:ibetahi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      REAL*8 alpha
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getdbgphipoint(dbgphi(i,j),
     &        beta(i,j),freq,x,alpha)
      enddo
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
      REAL*8 dbgphi
      REAL*8 beta
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 alpha
        REAL*8 gradphi(0:2 -1),gradbeta(0:2 -1)
        REAL*8 alphaphiplusbetalapphi,gradbetadotgradphi
        call getbetapoint(beta,freq,x)
        call getlofphipoint(alphaphiplusbetalapphi,freq,x,alpha,beta)
        call getgradbetapoint(gradbeta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        gradbetadotgradphi = gradbeta(0)*gradphi(0)
     &                               + gradbeta(1)*gradphi(1)
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
      REAL*8 beta
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
        beta = x(0)*x(0)
     &                 + x(1)*x(1)
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
      REAL*8 gradbeta(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
        integer idir
        do idir = 0, 2 -1
            gradbeta(idir) = (2.0d0)*x(idir)
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
      REAL*8 gradphi(0:1)
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
        integer idir
        REAL*8 beta
        call getbetapoint(beta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        do idir = 0, 2 -1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETSRC(
     &           src
     &           ,isrclo0,isrclo1
     &           ,isrchi0,isrchi1
     &           ,freq
     &           ,dx
     &           ,diffconst
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer isrclo0,isrclo1
      integer isrchi0,isrchi1
      REAL*8 src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1)
      REAL*8 freq(0:1)
      REAL*8 dx(0:1)
      REAL*8 diffconst
      REAL*8 problo(0:1)
      REAL*8 probhi(0:1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getsrcpoint(src(i,j),freq,x,diffconst)
      enddo
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
      REAL*8 src
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 diffconst
        REAL*8 fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2)
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1)))
        src = (-fac*diffconst)*phi
        return
        end