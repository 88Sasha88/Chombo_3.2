      subroutine BASEFABINTPLUS(
     &           sum
     &           ,isumlo0
     &           ,isumhi0
     &           ,nsumcomp
     &           ,piece
     &           ,ipiecelo0
     &           ,ipiecehi0
     &           ,npiececomp
     &           ,ibxlo0
     &           ,ibxhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nsumcomp
      integer isumlo0
      integer isumhi0
      integer sum(
     &           isumlo0:isumhi0,
     &           0:nsumcomp-1)
      integer npiececomp
      integer ipiecelo0
      integer ipiecehi0
      integer piece(
     &           ipiecelo0:ipiecehi0,
     &           0:npiececomp-1)
      integer ibxlo0
      integer ibxhi0
      integer i0
      integer var
      do var = 0, nsumcomp-1
      do i0 = ibxlo0,ibxhi0
            sum(i0, var) = sum(i0, var) +
     &        piece(i0, var)
      enddo
      enddo
      return
      end
