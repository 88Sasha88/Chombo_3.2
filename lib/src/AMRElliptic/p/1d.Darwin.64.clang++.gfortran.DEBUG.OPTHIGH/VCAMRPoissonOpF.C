#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
#if CH_SPACEDIM == 1
      subroutine GSRBHELMHOLTZVC1D(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,redBlack
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer redBlack
#elif CH_SPACEDIM == 2
      subroutine GSRBHELMHOLTZVC2D(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,redBlack
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer redBlack
#elif CH_SPACEDIM == 3
      subroutine GSRBHELMHOLTZVC3D(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,bCoef2
     &           ,ibCoef2lo0
     &           ,ibCoef2hi0
     &           ,nbCoef2comp
     &           ,lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,redBlack
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0
      integer ibCoef2hi0
      REAL_T bCoef2(
     &           ibCoef2lo0:ibCoef2hi0,
     &           0:nbCoef2comp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer redBlack
#else
      Somthing_that_will_not_compile
#endif
      REAL_T dxinv,lofphi
      integer n,ncomp,indtot,imin,imax
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = one/(dx*dx)
      do n = 0, ncomp - 1
#if CH_SPACEDIM==3
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              lofphi =
     &            alpha * aCoef(i,n) * phi(i,n)
     &          - beta  *
     &             (
     &               bCoef0(i+1,n)
     &               * (phi(i+1,n)-phi(i  ,n))
     &
     &             - bCoef0(i  ,n)
     &               * (phi(i  ,n)-phi(i-1,n)) 
     &             ) * dxinv
              phi(i,n) = phi(i,n)
     &          - lambda(i,n) * (lofphi - rhs(i,n))
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM==3
        enddo
#endif
      enddo
      return
      end
#if CH_SPACEDIM == 1
      subroutine VCCOMPUTEOP1D(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0
      integer ilofphihi0
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 2
      subroutine VCCOMPUTEOP2D(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0
      integer ilofphihi0
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 3
      subroutine VCCOMPUTEOP3D(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,bCoef2
     &           ,ibCoef2lo0
     &           ,ibCoef2hi0
     &           ,nbCoef2comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlofphicomp
      integer ilofphilo0
      integer ilofphihi0
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0
      integer ibCoef2hi0
      REAL_T bCoef2(
     &           ibCoef2lo0:ibCoef2hi0,
     &           0:nbCoef2comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#else
      Somthing_that_will_not_compile
#endif
      REAL_T dxinv
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          lofphi(i,n) =
     &        alpha * aCoef(i,n) * phi(i,n)
     &      - beta  *
     &         (
     &           bCoef0(i+1,n)
     &           * (phi(i+1,n) - phi(i  ,n))
     &
     &         - bCoef0(i  ,n)
     &           * (phi(i  ,n) - phi(i-1,n)) 
     &         ) * dxinv
        
      enddo
      enddo
      return
      end
#if CH_SPACEDIM == 1
      subroutine VCCOMPUTERES1D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 2
      subroutine VCCOMPUTERES2D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 3
      subroutine VCCOMPUTERES3D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,bCoef2
     &           ,ibCoef2lo0
     &           ,ibCoef2hi0
     &           ,nbCoef2comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0
      integer ibCoef2hi0
      REAL_T bCoef2(
     &           ibCoef2lo0:ibCoef2hi0,
     &           0:nbCoef2comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#else
      Somthing_that_will_not_compile
#endif
      REAL_T dxinv
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nrescomp) then
         call MAYDAYERROR()
      endif
      
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          res(i,n) =
     &        rhs(i,n)
     &      - (alpha * aCoef(i,n) * phi(i,n)
     &       - beta  *
     &          (
     &            bCoef0(i+1,n)
     &            * (phi(i+1,n) - phi(i  ,n))
     &
     &          - bCoef0(i  ,n)
     &            * (phi(i  ,n) - phi(i-1,n)) 
     &          ) * dxinv
     &        )
        
      enddo
      enddo
      return
      end
#if CH_SPACEDIM == 1
      subroutine RESTRICTRESVC1D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 2
      subroutine RESTRICTRESVC2D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#elif CH_SPACEDIM == 3
      subroutine RESTRICTRESVC3D(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,alpha
     &           ,aCoef
     &           ,iaCoeflo0
     &           ,iaCoefhi0
     &           ,naCoefcomp
     &           ,beta
     &           ,bCoef0
     &           ,ibCoef0lo0
     &           ,ibCoef0hi0
     &           ,nbCoef0comp
     &           ,bCoef1
     &           ,ibCoef1lo0
     &           ,ibCoef1hi0
     &           ,nbCoef1comp
     &           ,bCoef2
     &           ,ibCoef2lo0
     &           ,ibCoef2hi0
     &           ,nbCoef2comp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL_T alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL_T aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL_T beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL_T bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0
      integer ibCoef1hi0
      REAL_T bCoef1(
     &           ibCoef1lo0:ibCoef1hi0,
     &           0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0
      integer ibCoef2hi0
      REAL_T bCoef2(
     &           ibCoef2lo0:ibCoef2hi0,
     &           0:nbCoef2comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
#else
      Somthing_that_will_not_compile
#endif
      REAL_T denom,dxinv,lofphi
      integer n,ncomp
      integer i
      integer ii
      ncomp = nphicomp
      dxinv = one / (dx*dx)
      denom = D_TERM(2, *2, *2)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          
          ii = i/2 
          lofphi =
     &        alpha * aCoef(i,n) * phi(i,n)
     &      - beta  *
     &         (
     &           bCoef0(i+1,n)
     &           * (phi(i+1,n)-phi(i  ,n))
     &
     &         - bCoef0(i  ,n)
     &           * (phi(i  ,n)-phi(i-1,n)) 
     &         ) * dxinv
          res(ii,n) = res(ii,n)
     &                            + (rhs(i,n) - lofphi) / denom
        
      enddo
      enddo
      return
      end
      subroutine SUMFACES(
     &           lhs
     &           ,ilhslo0
     &           ,ilhshi0
     &           ,nlhscomp
     &           ,beta
     &           ,bCoefs
     &           ,ibCoefslo0
     &           ,ibCoefshi0
     &           ,nbCoefscomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,dir
     &           ,scale
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlhscomp
      integer ilhslo0
      integer ilhshi0
      REAL_T lhs(
     &           ilhslo0:ilhshi0,
     &           0:nlhscomp-1)
      REAL_T beta
      integer nbCoefscomp
      integer ibCoefslo0
      integer ibCoefshi0
      REAL_T bCoefs(
     &           ibCoefslo0:ibCoefshi0,
     &           0:nbCoefscomp-1)
      integer iboxlo0
      integer iboxhi0
      integer dir
      REAL_T scale
      REAL_T sumVal
      integer i
      integer ii
      integer n
      
      ii = CHF_ID(0,dir)
      do n = 0, nlhscomp-1
        
      do i = iboxlo0,iboxhi0

          sumVal = bCoefs(i+ii,n)
     &           + bCoefs(i   ,n)
          lhs(i,n) = lhs(i,n) + scale * beta * sumVal
        
      enddo
      enddo
      return
      end
