      subroutine AVERAGE(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,refRatio
     &           ,ibreflo0
     &           ,ibrefhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL*8 fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer iboxlo0
      integer iboxhi0
      integer refRatio
      integer ibreflo0
      integer ibrefhi0
      integer var
      integer ic0
      integer ip0
      integer ii0
      REAL*8 refScale,coarseSum
      refScale = (1.0d0) / (refRatio**1)
      if (refRatio .eq. 2) then
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = 2*ic0
            coarse(ic0,var) = refScale *
     &        (
     &          fine(ip0,var)
     &        + fine(ip0+1,var))
      enddo
        enddo
      else if (refRatio .eq. 4) then
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = 4*ic0
            coarse(ic0,var) = refScale *
     &        (
     &           fine(ip0  ,var)
     &         + fine(ip0+1,var)
     &         + fine(ip0+2,var)
     &         + fine(ip0+3,var) )
      enddo
        enddo
      else
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = ic0*refRatio
            coarseSum = 0
      do ii0 = ibreflo0,ibrefhi0
              coarseSum = coarseSum + fine( ip0+ii0,var)
      enddo
            coarse(ic0,var) = coarseSum * refScale
      enddo
       enddo
      endif
      return
      end
      subroutine AVERAGEHARMONIC(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,refRatio
     &           ,ibreflo0
     &           ,ibrefhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL*8 fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer iboxlo0
      integer iboxhi0
      integer refRatio
      integer ibreflo0
      integer ibrefhi0
      integer var
      integer ic0
      integer ip0
      integer ii0
      REAL*8 refScale,coarseSum
      refScale = (1.0d0) / (refRatio**1)
      if (refRatio .eq. 2) then
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = 2*ic0
            coarse(ic0,var) = refScale *
     &        (
     &          (1.0d0)/fine(ip0  ,var)
     &        + (1.0d0)/fine(ip0+1,var))
            coarse(ic0,var) = (1.0d0)/coarse(ic0,var)
      enddo
        enddo
      else if (refRatio .eq. 4) then
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = 4*ic0
            coarse(ic0,var) = refScale *
     &        (
     &           (1.0d0)/fine(ip0  ,var)
     &         + (1.0d0)/fine(ip0+1,var)
     &         + (1.0d0)/fine(ip0+2,var)
     &         + (1.0d0)/fine(ip0+3,var) )
            coarse(ic0,var) = (1.0d0)/coarse(ic0,var)
      enddo
        enddo
      else
        do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
            ip0 = ic0*refRatio
            coarseSum = 0
      do ii0 = ibreflo0,ibrefhi0
              coarseSum = coarseSum + (1.0d0)/fine( ip0+ii0,var)
      enddo
            coarse(ic0,var) = (1.0d0)/(coarseSum * refScale)
      enddo
       enddo
      endif
      return
      end
      subroutine AVERAGEINTVECTREF(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,ref
     &           ,ibreflo0
     &           ,ibrefhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL*8 fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer iboxlo0
      integer iboxhi0
      integer ref(0:0)
      integer ibreflo0
      integer ibrefhi0
      integer var
      integer ic0
      integer ip0
      integer ii0
      REAL*8 refScale,coarseSum
      refScale = (1.0d0) / (
     &  ref(0))
      do var = 0, ncoarsecomp - 1
      do ic0 = iboxlo0,iboxhi0
          ip0 = ic0*ref(0)
          coarseSum = 0
      do ii0 = ibreflo0,ibrefhi0
            coarseSum = coarseSum + fine( ip0+ii0,var)
      enddo
          coarse(ic0,var) = coarseSum * refScale
      enddo
      enddo
      return
      end
