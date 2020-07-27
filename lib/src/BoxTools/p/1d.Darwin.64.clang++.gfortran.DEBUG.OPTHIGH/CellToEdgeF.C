#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine CELLTOEDGE(
     &           cellData
     &           ,icellDatalo0
     &           ,icellDatahi0
     &           ,edgeData
     &           ,iedgeDatalo0
     &           ,iedgeDatahi0
     &           ,iedgeBoxlo0
     &           ,iedgeBoxhi0
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer icellDatalo0
      integer icellDatahi0
      REAL_T cellData(
     &           icellDatalo0:icellDatahi0)
      integer iedgeDatalo0
      integer iedgeDatahi0
      REAL_T edgeData(
     &           iedgeDatalo0:iedgeDatahi0)
      integer iedgeBoxlo0
      integer iedgeBoxhi0
      integer dir
      integer i
      integer ii
      
      do i = iedgeBoxlo0,iedgeBoxhi0

      
        ii = i-CHF_ID(0,dir)
        edgeData(i) = half*(
     &           cellData(ii)
     &         + cellData(i) )
        
      enddo
        return
        end
