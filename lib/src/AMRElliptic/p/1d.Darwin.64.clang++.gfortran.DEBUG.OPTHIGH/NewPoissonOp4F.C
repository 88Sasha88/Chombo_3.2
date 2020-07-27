#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine OPERATORLAP4(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           ,alpha
     &           ,beta
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
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      REAL_T dxinv, lap
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,n) - phi(i-2,n)
     & + sixteen*phi(i+1,n) - phi(i+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,n) )
     &       * twelfth * dxinv
          lofphi(i,n) = alpha*phi(i,n)+beta*lap
        
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES4(
     &           r
     &           ,irlo0
     &           ,irhi0
     &           ,nrcomp
     &           ,phi
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
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrcomp
      integer irlo0
      integer irhi0
      REAL_T r(
     &           irlo0:irhi0,
     &           0:nrcomp-1)
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
      REAL_T beta
      REAL_T dxinv, lap, lhs
      integer n,ncomp
      integer i
      ncomp = nphicomp
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,n) - phi(i-2,n)
     & + sixteen*phi(i+1,n) - phi(i+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,n) )
     &       * twelfth * dxinv
          lhs = alpha*phi(i,n) + beta*lap
          r(i,n) = rhs(i,n) - lhs
         
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES4(
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
     &           ,beta
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
      REAL_T beta
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T denom,dxinv,lofphi,lap
      integer n,ncomp
      integer i
      integer ii
      ncomp = nphicomp
      dxinv = one / (dx*dx)
      denom = D_TERM(2, *2, *2)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          
          ii = i/2 
          lap = ( 
     &   sixteen*phi(i-1,n) - phi(i-2,n)
     & + sixteen*phi(i+1,n) - phi(i+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,n) )
     &       * twelfth * dxinv
          lofphi = alpha*phi(i,n) + beta*lap
          res(ii,n) = res(ii,n)
     &                            + (rhs(i,n) - lofphi) / denom
        
      enddo
      enddo
      return
      end
      subroutine GSRBLAPLACIAN4(
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
     &           ,tmp
     &           ,itmplo0
     &           ,itmphi0
     &           ,ntmpcomp
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
      integer ntmpcomp
      integer itmplo0
      integer itmphi0
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           0:ntmpcomp-1)
      integer redBlack
      REAL_T dx2t, thD
      integer i
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               tmp(i,n) = thD*( 
     &           sixteen*phi(i+1,n) - phi(i+2,n)
     &         + sixteen*phi(i-1,n) - phi(i-2,n)
     &         - dx2t*rhs(i,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = thD*( 
     &           sixteen*tmp(i+1,n) - tmp(i+2,n)
     &         + sixteen*tmp(i-1,n) - tmp(i-2,n)
     &         - dx2t*rhs(i,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = tmp(i,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine SORLAPLACIAN4(
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
      REAL_T dx2t, thD, tmp, omega
      integer i
      integer n,ncomp
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      omega = 0.47
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         
      do i = iregionlo0,iregionhi0

         tmp = thD*( 
     &        sixteen*phi(i+1,n) - phi(i+2,n)
     &        + sixteen*phi(i-1,n) - phi(i-2,n)
     &        - dx2t*rhs(i,n) )
         phi(i,n) = omega*tmp
     &        + (one-omega)*phi(i,n)
         
      enddo
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ4(
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
     &           ,tmp
     &           ,itmplo0
     &           ,itmphi0
     &           ,ntmpcomp
     &           ,alpha
     &           ,beta
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
      integer ntmpcomp
      integer itmplo0
      integer itmphi0
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           0:ntmpcomp-1)
      REAL_T alpha
      REAL_T beta
      integer redBlack
      REAL_T dx2t, lambda, lap, dxinv, helm
      integer i
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      dxinv = one/(dx*dx)
      lambda = one/(alpha - beta*thirty*CH_SPACEDIM*twelfth*dxinv)
      lambda = lambda*(0.60)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
          lap = ( 
     &   sixteen*phi(i-1,n) - phi(i-2,n)
     & + sixteen*phi(i+1,n) - phi(i+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,n) )
     &       * twelfth * dxinv
          helm = alpha*phi(i,n) + beta*lap
          tmp(i,n) = phi(i,n) +
     &      lambda*( rhs(i,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               lap = ( 
     &           sixteen*tmp(i+1,n) - tmp(i+2,n)
     &         + sixteen*tmp(i-1,n) - tmp(i-2,n)
     &                     -(thirty*CH_SPACEDIM)*tmp(i,n) )
     &       * twelfth * dxinv
               helm = alpha*tmp(i,n) + beta*lap
               phi(i,n) = tmp(i,n) +
     &              lambda*( rhs(i,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = tmp(i,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine NEWGETFLUX4(
     &           flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,beta_dx
     &           ,a_idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iboxlo0
      integer iboxhi0
      REAL_T beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii
      integer i 
      ncomp = nphicomp
      
      ii = CHF_ID(a_idir, 0)
      do n = 0, ncomp-1
          
      do i = iboxlo0,iboxhi0

          flux(i,n) = beta_dx * twelfth *
     &        ( fifteen*phi(i,n)
     &           + phi(i-2*ii,n)
     &           - phi(i+ii,n)
     &           - fifteen*phi(i-ii,n) )
          
      enddo
      enddo
      return
      end
      subroutine PROLONGLINEAR(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,ifineBoxlo0
     &           ,ifineBoxhi0
     &           ,icrseBoxlo0
     &           ,icrseBoxhi0
     &           ,r
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
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer ifineBoxlo0
      integer ifineBoxhi0
      integer icrseBoxlo0
      integer icrseBoxhi0
      integer r
      INTEGER ncomp, n
      integer i 
      integer ic
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do i = ifineBoxlo0,ifineBoxhi0

          
           ic = i/r
           phi(i,n) =  phi(i,n) +
     &          coarse(ic,n)
          
           if (ic.ne.icrseBoxhi0 .and.
     &         (ic*r.lt.i .or. ic.eq.icrseBoxlo0)) then
              phi(i,n) =  phi(i,n) +
     &             (coarse(ic+1,n)
     &              - coarse(ic,n))/r*(i+half-ic*r-half*r)
           else
              phi(i,n) =  phi(i,n) +
     &             (- coarse(ic-1,n)
     &              + coarse(ic,n))/r*(i+half-ic*r-half*r)
           endif
          
      enddo
      enddo
      return
      end
