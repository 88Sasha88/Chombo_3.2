#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEOPLAP(
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
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      do var = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

            
            lphi = ( (phi(i+1 , var)
     &              - phi(i   , var) )
     &            -  (phi(i   , var)
     &              - phi(i-1 , var) ) ) * dxinv2 
            lofphi(i, var) =  lphi
         
      enddo
      end do
      return
      end
      subroutine NODEOPLAPPOINT(
     &           lofphi
     &           ,ilofphilo0
     &           ,ilofphihi0
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,pt
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
      integer pt(0:0)
      REAL_T dx
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      
      i = pt(0) 
      do var = 0, ncomp-1
         
         lphi = ( (phi(i+1 , var)
     &           - phi(i   , var) )
     &         -  (phi(i   , var)
     &           - phi(i-1 , var) ) ) * dxinv2 
         lofphi(i, var) =  lphi
      end do
      return
      end
      subroutine NODEGRAD(
     &           grdphi
     &           ,igrdphilo0
     &           ,igrdphihi0
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0
      integer igrdphihi0
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
         
      do i = iregionlo0,iregionhi0

            
            grdphi(i, gbase) =
     &           ( phi(i+1 , var)
     &           - phi(i-1 , var) ) * dxinvh 
         
      enddo
      end do
      return
      end
      subroutine NODEGRADPOINT(
     &           grdphi
     &           ,igrdphilo0
     &           ,igrdphihi0
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,pt
     &           ,dx
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ngrdphicomp
      integer igrdphilo0
      integer igrdphihi0
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer pt(0:0)
      REAL_T dx
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      
      i = pt(0) 
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
            
            grdphi(i, gbase) =
     &           ( phi(i+1 , var)
     &           - phi(i-1 , var) ) * dxinvh 
      end do
      return
      end
      subroutine NODEGSRBLEVELLAP(
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
      REAL_T lambda
      REAL_T dxinv2, lphi
      integer i
      integer imin, imax, var, ncomp, indtot
      dxinv2 = one/(dx*dx)
      lambda = (dx*dx) / (two*CH_SPACEDIM)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAY_ERROR()
      endif
      do var = 0, ncomp - 1
#if CH_SPACEDIM>=3
         do k = iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM>=2
            do j = iregionlo1, iregionhi1
#endif
               imin = iregionlo0
               indtot = imin 
               imin = imin + mod(indtot + redBlack, 2)
               imax = iregionhi0
               do i = imin, imax, 2
#ifndef NDEBUG
                  if (mod(i, 2) .ne. redBlack) then
                     print *, 'NODEGSRBLEVELLAP:  computing ',
     &                    i, 
     &                    ' at pass ', redBlack
                  endif
#endif
            
            lphi = ( (phi(i+1 , var)
     &              + phi(i-1 , var) ) ) * dxinv2 
                  phi(i, var) =
     &                 lambda * (lphi - rhs(i, var))
               enddo
#if CH_SPACEDIM>=2
            enddo
#endif
#if CH_SPACEDIM>=3
         enddo
#endif
      enddo
      return
      end
