#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
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
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0)
      REAL_T freq(0:0)
      REAL_T dx(0:0)
      REAL_T problo(0:0)
      REAL_T probhi(0:0)
      REAL_T aCoef
      REAL_T bCoef
      integer iboxlo0
      integer iboxhi0
        integer i
        real_t x(0:CH_SPACEDIM-1)
        
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
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
      REAL_T lofphi
      REAL_T freq(0:0)
      REAL_T x(0:0)
      REAL_T aCoefmult
      REAL_T bCoefmult
        integer dir
        real_t fac,phi, temp
        fac = -(freq(0)**2)
        fac = fac *(x(0))
        phi = (sin(freq(0)*x(0)))
        lofphi = fac*phi
        temp = 0.0
        do dir=0, CH_SPACEDIM-1
           if (dir.eq.0) then
              temp = freq(0)*cos(freq(0)*x(0))
#if CH_SPACEDIM > 1
           else if (dir.eq.1) then
              temp = freq(1)*sin(freq(0)*x(0))
#if CH_SPACEDIM > 2
           else if (dir.eq.2) then
              temp = freq(2)*sin(freq(0)*x(0))
#endif
#endif
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
      REAL_T gradphi(0:0)
      REAL_T freq(0:0)
      REAL_T x(0:0)
        
        gradphi(0) = freq(0) * cos(freq(0)*x(0))                                       
        return
        end
