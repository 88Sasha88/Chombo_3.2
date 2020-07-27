#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine DOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0
     &           ,iafabhi0
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0
     &           ,ibfabhi0
     &           ,nbfabcomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0
      integer iafabhi0
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0
      integer ibfabhi0
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           0:nbfabcomp-1)
      integer iregionlo0
      integer iregionhi0
      integer startcomp
      integer endcomp
      integer i0
      integer nv
      dotprodout = 0
      do nv=startcomp,endcomp,1
         
      do i0 = iregionlo0,iregionhi0

         dotprodout = dotprodout +
     &        afab(i0,nv)*
     &        bfab(i0,nv)
         
      enddo
      enddo
      return
      end
