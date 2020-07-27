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
      REAL*8 dxinv2, lphi
      integer var, ncomp
      integer  i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = (1.0d0) / (dx*dx)
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
      REAL*8 lofphi(
     &           ilofphilo0:ilofphihi0,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer pt(0:0)
      REAL*8 dx
      REAL*8 dxinv2, lphi
      integer var, ncomp
      integer  i
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = (1.0d0) / (dx*dx)
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
      REAL*8 grdphi(
     &           igrdphilo0:igrdphihi0,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 dxinvh
      integer var, ncomp, gbase
      integer  i
      ncomp = nphicomp
      if (1 * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = (0.500d0) / dx
      do var = 0, ncomp-1
         gbase = 1 * var
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
      REAL*8 grdphi(
     &           igrdphilo0:igrdphihi0,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer pt(0:0)
      REAL*8 dx
      REAL*8 dxinvh
      integer var, ncomp, gbase
      integer  i
      ncomp = nphicomp
      if (1 * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = (0.500d0) / dx
      i = pt(0) 
      do var = 0, ncomp-1
         gbase = 1 * var
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
      integer redBlack
      REAL*8 lambda
      REAL*8 dxinv2, lphi
      integer i
      integer imin, imax, var, ncomp, indtot
      dxinv2 = (1.0d0)/(dx*dx)
      lambda = (dx*dx) / ((2.0d0)*1)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAY_ERROR()
      endif
      do var = 0, ncomp - 1
               imin = iregionlo0
               indtot = imin 
               imin = imin + mod(indtot + redBlack, 2)
               imax = iregionhi0
               do i = imin, imax, 2
            lphi = ( (phi(i+1 , var)
     &              + phi(i-1 , var) ) ) * dxinv2 
                  phi(i, var) =
     &                 lambda * (lphi - rhs(i, var))
               enddo
      enddo
      return
      end
