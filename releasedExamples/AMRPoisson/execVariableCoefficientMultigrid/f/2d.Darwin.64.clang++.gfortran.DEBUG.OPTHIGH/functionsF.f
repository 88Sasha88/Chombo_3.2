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
        subroutine GETLOFPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1
     &           ,ilofphihi0,ilofphihi1
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,aCoef
     &           ,bCoef
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
      REAL*8 aCoef
      REAL*8 bCoef
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
        integer i,j
        REAL*8 x(0:2 -1)
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          x(1) = (j+(0.500d0))*dx(1) + problo(1)
          call getlofphipoint(lofphi(i,j),freq,x,aCoef,bCoef)
      enddo
      enddo
        return
        end
        subroutine GETLOFPHIPOINT(
     &           lofphi
     &           ,freq
     &           ,x
     &           ,aCoefmult
     &           ,bCoefmult
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 lofphi
      REAL*8 freq(0:1)
      REAL*8 x(0:1)
      REAL*8 aCoefmult
      REAL*8 bCoefmult
        integer dir
        REAL*8 fac,phi, temp
        fac = -(freq(0)**2
     &                  + freq(1)**2)
        fac = fac *(x(0) +x(1))
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1)))
        lofphi = fac*phi
        temp = 0.0
        do dir=0, 2 -1
           if (dir.eq.0) then
              temp = freq(0)*cos(freq(0)*x(0))
     &                                *sin(freq(1)*x(1))
           else if (dir.eq.1) then
              temp = freq(1)*sin(freq(0)*x(0))
     &                                *cos(freq(1)*x(1))
           endif
           lofphi = lofphi + temp
        enddo
        lofphi = aCoefmult*x(0)*phi + bCoefmult*lofphi
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
