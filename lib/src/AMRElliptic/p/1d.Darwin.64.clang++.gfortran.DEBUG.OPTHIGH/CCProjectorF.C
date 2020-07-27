#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine CCPEXTRAPTODOMFACE(
     &           facevel
     &           ,ifacevello0
     &           ,ifacevelhi0
     &           ,facedir
     &           ,ishift
     &           ,ifaceboxlo0
     &           ,ifaceboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifacevello0
      integer ifacevelhi0
      REAL_T facevel(
     &           ifacevello0:ifacevelhi0)
      integer facedir
      integer ishift
      integer ifaceboxlo0
      integer ifaceboxhi0
      integer i
      integer ioff
      
      ioff = ishift*chf_id(0,facedir)
      
      do i = ifaceboxlo0,ifaceboxhi0

      facevel(i) =
     &     (two*facevel(i+  ioff)
     &     -    facevel(i+2*ioff))
      
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1) 
      return
      end
      subroutine CCPAVECELLTOFACE(
     &           facevel
     &           ,ifacevello0
     &           ,ifacevelhi0
     &           ,cellvel
     &           ,icellvello0
     &           ,icellvelhi0
     &           ,facedir
     &           ,ifaceboxlo0
     &           ,ifaceboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifacevello0
      integer ifacevelhi0
      REAL_T facevel(
     &           ifacevello0:ifacevelhi0)
      integer icellvello0
      integer icellvelhi0
      REAL_T cellvel(
     &           icellvello0:icellvelhi0)
      integer facedir
      integer ifaceboxlo0
      integer ifaceboxhi0
      integer i
      integer ioff
      
      ioff = chf_id(0,facedir)
      
      do i = ifaceboxlo0,ifaceboxhi0

      facevel(i) =
     &     ( cellvel(i     )
     &     + cellvel(i-ioff)
     &     )*half
      
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*2  
      return
      end
      subroutine CCPAVEFACETOCELL(
     &           cellgrad
     &           ,icellgradlo0
     &           ,icellgradhi0
     &           ,facegrad
     &           ,ifacegradlo0
     &           ,ifacegradhi0
     &           ,facedir
     &           ,icellboxlo0
     &           ,icellboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icellgradlo0
      integer icellgradhi0
      REAL_T cellgrad(
     &           icellgradlo0:icellgradhi0)
      integer ifacegradlo0
      integer ifacegradhi0
      REAL_T facegrad(
     &           ifacegradlo0:ifacegradhi0)
      integer facedir
      integer icellboxlo0
      integer icellboxhi0
      integer i
      integer ioff
      
      ioff = chf_id(0,facedir)
      
      do i = icellboxlo0,icellboxhi0

      cellgrad(i) =
     &     (facegrad(i+ioff)
     &     +facegrad(i     ))*half
      
      enddo
      ch_flops=ch_flops+(icellboxhi0- icellboxlo0+1)*2  
      return
      end
      subroutine CCPCELLGRADFROMFACEDATA(
     &           cellgrad
     &           ,icellgradlo0
     &           ,icellgradhi0
     &           ,facedata
     &           ,ifacedatalo0
     &           ,ifacedatahi0
     &           ,facedir
     &           ,dx
     &           ,icellboxlo0
     &           ,icellboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icellgradlo0
      integer icellgradhi0
      REAL_T cellgrad(
     &           icellgradlo0:icellgradhi0)
      integer ifacedatalo0
      integer ifacedatahi0
      REAL_T facedata(
     &           ifacedatalo0:ifacedatahi0)
      integer facedir
      REAL_T dx
      integer icellboxlo0
      integer icellboxhi0
      integer i
      integer ioff
      
      ioff = chf_id(0,facedir)
      
      do i = icellboxlo0,icellboxhi0

      cellgrad(i) =
     &     (facedata(i+ioff)
     &     -facedata(i     ))/dx
      
      enddo
      ch_flops=ch_flops+(icellboxhi0- icellboxlo0+1)*2  
      return
      end
