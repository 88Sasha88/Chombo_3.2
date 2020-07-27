#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine GSMCAMRPOP(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,icoloredboxlo0
     &           ,icoloredboxhi0
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0)
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      real_t lambda, dxinv, sum_b, lphi
      integer i
      integer idir
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = one/(alpha - beta*sum_b)
      
      do i = icoloredBoxlo0,icoloredBoxhi0,2

        lphi = 
     &     (    phi(i+1)
     &     +    phi(i-1)
     $     -two*phi(i  )) 
        lphi = lphi*dxinv
        phi(i) =
     $       phi(   i) +
     &       lambda*(   rhs(   i) - lphi)
      
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ(
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
      REAL_T alpha
      REAL_T beta
      integer redBlack
      REAL_T lambda, dxinv, sum_b, lphi, helmop
      integer i
      integer n,ncomp,idir,indtot,imin,imax
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/(alpha - beta*sum_b)
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
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
              
              lphi =   (phi(i+1,n)
     &             +    phi(i-1,n)
     &             -two*CH_SPACEDIM*phi(i,n))*dxinv
              helmop = alpha*phi(i,n) + beta*lphi
              phi(i,n) = phi(i,n) +
     &             lambda*(helmop - rhs(i,n))
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
      subroutine GSRBLAPLACIAN(
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
      integer redBlack
      REAL_T lambda, dxinv, sum_b, lphi, lap
      integer i
      integer n,ncomp,idir,indtot,imin,imax,counter
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/sum_b
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
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
#if 0
              
              lphi =  ((phi(i+1,n)
     &                   - phi(i,n))
     &                   -(phi(i,n)
     &                   - phi(i-1,n)))*dxinv
#else
              lphi = ( 
     &                  phi(i+1,n)
     &             +    phi(i-1,n)
     &             ) * dxinv
#endif
              phi(i,n) = lambda*(rhs(i,n)-lphi)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM==3
        enddo
#endif
      enddo
      i=(iregionhi0- iregionlo0+1)
      ch_flops=ch_flops+ ncomp*i/2*(2*CH_SPACEDIM+3)
      return
      end
      subroutine GSRBLAZY(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,lphi
     &           ,ilphilo0
     &           ,ilphihi0
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,icoloredboxlo0
     &           ,icoloredboxhi0
     &           ,alpha
     &           ,beta
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0)
      integer ilphilo0
      integer ilphihi0
      REAL_T lphi(
     &           ilphilo0:ilphihi0)
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      REAL_T alpha
      REAL_T beta
      REAL_T dx
      integer i, idir
      real_t dxinv, sum_b, lambda
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/(alpha - beta*sum_b)
      
      do i = icoloredBoxlo0,icoloredBoxhi0,2

      phi(i) =
     $     phi(   i) -
     &     lambda*(
     $     rhs(   i) -
     $     lphi(  i))
      
      enddo
      return
      end
      subroutine AMRPMULTICOLOR(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,weight
     &           ,alpha
     &           ,beta
     &           ,dx
     &           ,icoloredboxlo0
     &           ,icoloredboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0)
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0)
      REAL_T weight
      REAL_T alpha
      REAL_T beta
      REAL_T dx(0:0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      integer i
      real_t laplphi, dx0
      dx0 = beta/(dx(0) * dx(0))
      
      do i = icoloredBoxlo0,icoloredBoxhi0,2

        laplphi = 
     &     (    phi(i+1)
     &     +    phi(i-1)
     $     -two*phi(i  ))*dx0 
        laplphi = laplphi + alpha * phi(i)
        phi(i) = phi(i) +
     &     weight*(rhs(i) - laplphi)
      
      enddo
      return
      end
      subroutine OPERATORLAP(
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
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do i = iregionlo0,iregionhi0

          lap = (  phi(i+1,n)
     &                     + phi(i-1,n) 
     &                     -(2*CH_SPACEDIM)*phi(i,n) )
     &       * dxinv
          lofphi(i,n) = alpha*phi(i,n)+beta*lap
        
      enddo
      enddo
      i=(iregionhi0- iregionlo0+1)
      ch_flops=ch_flops+ ncomp*i*(2*CH_SPACEDIM+5)
      return
      end
      subroutine OPERATORLAPRES(
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
      REAL_T dxinv, lap
      integer n,ncomp
      integer i
      ncomp = nphicomp
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

          lap = (  phi(i+1,n)
     &                     + phi(i-1,n) 
     &                     -(2*CH_SPACEDIM)*phi(i,n) )
     &       * dxinv
         r(i,n) = -alpha*phi(i,n) -beta*lap +
     &       rhs(i,n)
         
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*ncomp*(2*CH_SPACEDIM+6)
      return
      end
      subroutine RESTRICTRES(
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
          lofphi = alpha * phi(i,n)
     &           + beta  *
     &              ( phi(i+1,n)
     &                        + phi(i-1,n) 
     &                        - phi(i  ,n) * 2 * CH_SPACEDIM
     &              ) * dxinv
          res(ii,n) = res(ii,n)
     &                            + (rhs(i,n) - lofphi) / denom
        
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*ncomp*(2*CH_SPACEDIM+6)
      return
      end
      subroutine PROLONG(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,m
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
      integer iregionlo0
      integer iregionhi0
      integer m
      INTEGER ncomp, n
      integer i
      integer ii
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do i = iregionlo0,iregionhi0

          
          ii = i/m
          phi(i,n) =  phi(i,n) +
     &        coarse(ii,n)
         
      enddo
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)
      return
      end
      subroutine PROLONG_2(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,m
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
      integer iregionlo0
      integer iregionhi0
      integer m
      INTEGER ncomp, n
      integer i
      integer offs(CH_SPACEDIM)
      integer ic
      real_t f0, den, fx(CH_SPACEDIM)
      den = one/(4**CH_SPACEDIM)
      
      fx(1) = three*den
      f0 = one*den
      ncomp = nphicomp
      
      do i = iregionlo0,iregionhi0

        
        ic = i/m
        
        offs(1) = 2*mod(i,2) - 1
        do n = 0, ncomp-1
          phi(i,n) = phi(i,n)
     $      + fx(CH_SPACEDIM)*
     $        coarse(ic,n)
     $      + f0*coarse(ic+offs(1),n) 
#if CH_SPACEDIM > 1
          phi(i,n) = phi(i,n) 
     $      + fx(CH_SPACEDIM-1)*
     $        (
     $         coarse(ic+offs(1),n)  )
#if CH_SPACEDIM > 2
          phi(i,n) = phi(i,n) 
     $      + fx(CH_SPACEDIM-2)*
     $        (
     $         coarse(ic+offs(1),n)  )
#endif
#endif
        enddo
      
      enddo
      return
      end
      subroutine NEWGETFLUX(
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

          flux(i,n) =
     &        (phi(i,n)-
     &         phi(i-ii,n))*beta_dx
          
      enddo
      enddo
      return
      end
