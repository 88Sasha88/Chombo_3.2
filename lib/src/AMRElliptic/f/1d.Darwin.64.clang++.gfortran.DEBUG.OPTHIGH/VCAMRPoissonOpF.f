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
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL*8 lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer redBlack
      REAL*8 dxinv,lofphi
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
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp - 1
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
     &             - bCoef0(i  ,n)
     &               * (phi(i  ,n)-phi(i-1,n)) 
     &             ) * dxinv
              phi(i,n) = phi(i,n)
     &          - lambda(i,n) * (lofphi - rhs(i,n))
            enddo
      enddo
      return
      end
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
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          lofphi(i,n) =
     &        alpha * aCoef(i,n) * phi(i,n)
     &      - beta  *
     &         (
     &           bCoef0(i+1,n)
     &           * (phi(i+1,n) - phi(i  ,n))
     &         - bCoef0(i  ,n)
     &           * (phi(i  ,n) - phi(i-1,n)) 
     &         ) * dxinv
      enddo
      enddo
      return
      end
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
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nrescomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif                                  
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          res(i,n) =
     &        rhs(i,n)
     &      - (alpha * aCoef(i,n) * phi(i,n)
     &       - beta  *
     &          (
     &            bCoef0(i+1,n)
     &            * (phi(i+1,n) - phi(i  ,n))
     &          - bCoef0(i  ,n)
     &            * (phi(i  ,n) - phi(i-1,n)) 
     &          ) * dxinv
     &        )
      enddo
      enddo
      return
      end
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
      REAL*8 res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0
      integer iaCoefhi0
      REAL*8 aCoef(
     &           iaCoeflo0:iaCoefhi0,
     &           0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0
      integer ibCoef0hi0
      REAL*8 bCoef0(
     &           ibCoef0lo0:ibCoef0hi0,
     &           0:nbCoef0comp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi
      integer n,ncomp
      integer i
      integer ii
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          ii = i/2 
          lofphi =
     &        alpha * aCoef(i,n) * phi(i,n)
     &      - beta  *
     &         (
     &           bCoef0(i+1,n)
     &           * (phi(i+1,n)-phi(i  ,n))
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
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlhscomp
      integer ilhslo0
      integer ilhshi0
      REAL*8 lhs(
     &           ilhslo0:ilhshi0,
     &           0:nlhscomp-1)
      REAL*8 beta
      integer nbCoefscomp
      integer ibCoefslo0
      integer ibCoefshi0
      REAL*8 bCoefs(
     &           ibCoefslo0:ibCoefshi0,
     &           0:nbCoefscomp-1)
      integer iboxlo0
      integer iboxhi0
      integer dir
      REAL*8 scale
      REAL*8 sumVal
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
