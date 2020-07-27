        subroutine GETPHIPOINT(
     &           phi
     &           ,freq
     &           ,x
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL*8 phi
      REAL*8 freq(0:0)
      REAL*8 x(0:0)
        phi = sin(freq(0)*x(0))
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
     &           ,aCoef
     &           ,bCoef
     &           ,iboxlo0
     &           ,iboxhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ilofphilo0
      integer ilofphihi0
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0)
      REAL*8 freq(0:0)
      REAL*8 dx(0:0)
      REAL*8 problo(0:0)
      REAL*8 probhi(0:0)
      REAL*8 aCoef
      REAL*8 bCoef
      integer iboxlo0
      integer iboxhi0
        integer i
        REAL*8 x(0:1 -1)
      do i = iboxlo0,iboxhi0
          x(0) = (i+(0.500d0))*dx(0) + problo(0)
          call getlofphipoint(lofphi(i),freq,x,aCoef,bCoef)
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
      REAL*8 freq(0:0)
      REAL*8 x(0:0)
      REAL*8 aCoefmult
      REAL*8 bCoefmult
        integer dir
        REAL*8 fac,phi, temp
        fac = -(freq(0)**2)
        fac = fac *(x(0))
        phi = (sin(freq(0)*x(0)))
        lofphi = fac*phi
        temp = 0.0
        do dir=0, 1 -1
           if (dir.eq.0) then
              temp = freq(0)*cos(freq(0)*x(0))
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
      REAL*8 gradphi(0:0)
      REAL*8 freq(0:0)
      REAL*8 x(0:0)
        gradphi(0) = freq(0) * cos(freq(0)*x(0))                        
        return
        end