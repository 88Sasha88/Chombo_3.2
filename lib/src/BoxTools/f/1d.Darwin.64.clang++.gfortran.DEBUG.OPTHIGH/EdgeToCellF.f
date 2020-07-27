      subroutine EDGETOCELL(
     &           edgeData
     &           ,iedgeDatalo0
     &           ,iedgeDatahi0
     &           ,cellData
     &           ,icellDatalo0
     &           ,icellDatahi0
     &           ,icellBoxlo0
     &           ,icellBoxhi0
     &           ,dir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iedgeDatalo0
      integer iedgeDatahi0
      REAL*8 edgeData(
     &           iedgeDatalo0:iedgeDatahi0)
      integer icellDatalo0
      integer icellDatahi0
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0)
      integer icellBoxlo0
      integer icellBoxhi0
      integer dir
      integer i
      integer ii
      do i = icellBoxlo0,icellBoxhi0
      ii = i+CHF_ID(0,dir)
      cellData(i) = (0.500d0)*(
     &                    edgeData(i)
     &                   +edgeData(ii))
      enddo
      return
      end
      subroutine EDGETOINCREMENTCELL(
     &           edgeData
     &           ,iedgeDatalo0
     &           ,iedgeDatahi0
     &           ,cellData
     &           ,icellDatalo0
     &           ,icellDatahi0
     &           ,icellBoxlo0
     &           ,icellBoxhi0
     &           ,dir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iedgeDatalo0
      integer iedgeDatahi0
      REAL*8 edgeData(
     &           iedgeDatalo0:iedgeDatahi0)
      integer icellDatalo0
      integer icellDatahi0
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0)
      integer icellBoxlo0
      integer icellBoxhi0
      integer dir
      integer i0
      integer ii0
      ii0=CHF_ID(0, dir)
      do i0 = icellBoxlo0,icellBoxhi0
         cellData(i0) = cellData(i0) + (0.500d0)*(
     &      edgeData(i0) + edgeData(i0+ii0))
      enddo
      return
      end
      subroutine EDGETOCELLMAX(
     &           edgeData
     &           ,iedgeDatalo0
     &           ,iedgeDatahi0
     &           ,cellData
     &           ,icellDatalo0
     &           ,icellDatahi0
     &           ,icellBoxlo0
     &           ,icellBoxhi0
     &           ,dir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer iedgeDatalo0
      integer iedgeDatahi0
      REAL*8 edgeData(
     &           iedgeDatalo0:iedgeDatahi0)
      integer icellDatalo0
      integer icellDatahi0
      REAL*8 cellData(
     &           icellDatalo0:icellDatahi0)
      integer icellBoxlo0
      integer icellBoxhi0
      integer dir
      integer i
      integer ii
      do i = icellBoxlo0,icellBoxhi0
      ii = i+CHF_ID(0,dir)
      cellData(i) = max(
     &                    edgeData(i),
     &                    edgeData(ii))
      enddo
      return
      end
