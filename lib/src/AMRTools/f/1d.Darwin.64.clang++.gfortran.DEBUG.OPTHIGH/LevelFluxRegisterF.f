      subroutine INCREMENTFINE(
     &           fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,cFine
     &           ,icFinelo0
     &           ,icFinehi0
     &           ,ncFinecomp
     &           ,ifineBoxlo0
     &           ,ifineBoxhi0
     &           ,nRef
     &           ,scale
     &           ,srcStart
     &           ,destStart
     &           ,ncomp
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL*8 fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer ncFinecomp
      integer icFinelo0
      integer icFinehi0
      REAL*8 cFine(
     &           icFinelo0:icFinehi0,
     &           0:ncFinecomp-1)
      integer ifineBoxlo0
      integer ifineBoxhi0
      integer nRef(0:0)
      REAL*8 scale
      integer srcStart
      integer destStart
      integer ncomp
      integer i0
      integer ii0
      integer var, srcComp, destComp
      do var=0, ncomp-1
         srcComp = srcStart + var
         destComp = destStart + var
      do i0 = ifineBoxlo0,ifineBoxhi0
            ii0=i0/nRef(0)
            cFine(ii0,destComp) =
     &           cFine(ii0, destComp) +
     &           scale * fine(i0, srcComp)
      enddo
      enddo
      return
      end
