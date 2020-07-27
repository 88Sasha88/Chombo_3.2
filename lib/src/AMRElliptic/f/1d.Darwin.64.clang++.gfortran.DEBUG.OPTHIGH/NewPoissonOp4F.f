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
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dxinv, lap
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          lap = ( 
     &   (16.0d0)*phi(i-1,n) - phi(i-2,n)
     & + (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &                     -((30.0d0)*1)*phi(i,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
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
      REAL*8 r(
     &           irlo0:irhi0,
     &           0:nrcomp-1)
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
      REAL*8 beta
      REAL*8 dxinv, lap, lhs
      integer n,ncomp
      integer i
      ncomp = nphicomp
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          lap = ( 
     &   (16.0d0)*phi(i-1,n) - phi(i-2,n)
     & + (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &                     -((30.0d0)*1)*phi(i,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
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
      REAL*8 beta
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi,lap
      integer n,ncomp
      integer i
      integer ii
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          ii = i/2 
          lap = ( 
     &   (16.0d0)*phi(i-1,n) - phi(i-2,n)
     & + (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &                     -((30.0d0)*1)*phi(i,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
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
      integer ntmpcomp
      integer itmplo0
      integer itmphi0
      REAL*8 tmp(
     &           itmplo0:itmphi0,
     &           0:ntmpcomp-1)
      integer redBlack
      REAL*8 dx2t, thD
      integer i
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = (12.0d0)*dx*dx
      thD  = (1.000d0 / 30.000d0)/1
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               tmp(i,n) = thD*( 
     &           (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &         + (16.0d0)*phi(i-1,n) - phi(i-2,n)
     &         - dx2t*rhs(i,n) )
            enddo
        else if (redBlack .eq. black) then
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = thD*( 
     &           (16.0d0)*tmp(i+1,n) - tmp(i+2,n)
     &         + (16.0d0)*tmp(i-1,n) - tmp(i-2,n)
     &         - dx2t*rhs(i,n) )
            enddo
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = tmp(i,n)
            enddo
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
      REAL*8 dx2t, thD, tmp, omega
      integer i
      integer n,ncomp
      dx2t = (12.0d0)*dx*dx
      thD  = (1.000d0 / 30.000d0)/1
      omega = 0.47
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
      do i = iregionlo0,iregionhi0
         tmp = thD*( 
     &        (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &        + (16.0d0)*phi(i-1,n) - phi(i-2,n)
     &        - dx2t*rhs(i,n) )
         phi(i,n) = omega*tmp
     &        + ((1.0d0)-omega)*phi(i,n)
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
      integer ntmpcomp
      integer itmplo0
      integer itmphi0
      REAL*8 tmp(
     &           itmplo0:itmphi0,
     &           0:ntmpcomp-1)
      REAL*8 alpha
      REAL*8 beta
      integer redBlack
      REAL*8 dx2t, lambda, lap, dxinv, helm
      integer i
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = (12.0d0)*dx*dx
      dxinv = (1.0d0)/(dx*dx)
      lambda = (1.0d0)/(alpha - beta*(30.0d0)*1*(1.000d0 / 12.000d0)*dxi
     &nv)
      lambda = lambda*(0.60)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
          lap = ( 
     &   (16.0d0)*phi(i-1,n) - phi(i-2,n)
     & + (16.0d0)*phi(i+1,n) - phi(i+2,n)
     &                     -((30.0d0)*1)*phi(i,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
          helm = alpha*phi(i,n) + beta*lap
          tmp(i,n) = phi(i,n) +
     &      lambda*( rhs(i,n) - helm )
            enddo
        else if (redBlack .eq. black) then
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               lap = ( 
     &           (16.0d0)*tmp(i+1,n) - tmp(i+2,n)
     &         + (16.0d0)*tmp(i-1,n) - tmp(i-2,n)
     &                     -((30.0d0)*1)*tmp(i,n) )
     &       * (1.000d0 / 12.000d0) * dxinv
               helm = alpha*tmp(i,n) + beta*lap
               phi(i,n) = tmp(i,n) +
     &              lambda*( rhs(i,n) - helm )
            enddo
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,n) = tmp(i,n)
            enddo
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
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL*8 flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iboxlo0
      integer iboxhi0
      REAL*8 beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii
      integer i 
      ncomp = nphicomp
      ii = CHF_ID(a_idir, 0)
      do n = 0, ncomp-1
      do i = iboxlo0,iboxhi0
          flux(i,n) = beta_dx * (1.000d0 / 12.000d0) *
     &        ( (15.0d0)*phi(i,n)
     &           + phi(i-2*ii,n)
     &           - phi(i+ii,n)
     &           - (15.0d0)*phi(i-ii,n) )
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
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
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
     &              - coarse(ic,n))/r*(i+(0.500d0)-ic*r-(0.500d0)*r)
           else
              phi(i,n) =  phi(i,n) +
     &             (- coarse(ic-1,n)
     &              + coarse(ic,n))/r*(i+(0.500d0)-ic*r-(0.500d0)*r)
           endif
      enddo
      enddo
      return
      end
