#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine HALVEINT(
     &           arr
     &           ,iarrlo0
     &           ,iarrhi0
     &           ,narrcomp
     &           ,ibxlo0
     &           ,ibxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer narrcomp
      integer iarrlo0
      integer iarrhi0
      integer arr(
     &           iarrlo0:iarrhi0,
     &           0:narrcomp-1)
      integer ibxlo0
      integer ibxhi0
      integer i0
      integer var
      do var = 0, narrcomp-1
         
      do i0 = ibxlo0,ibxhi0

            arr(i0, var) = arr(i0, var) / 2
         
      enddo
      enddo
      return
      end
