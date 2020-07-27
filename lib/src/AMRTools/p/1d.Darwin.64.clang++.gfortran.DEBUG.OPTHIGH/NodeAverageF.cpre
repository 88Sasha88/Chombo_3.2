      subroutine NODEAVERAGE(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0
     &           ,iweighthi0
     &           ,nweightcomp
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
      integer iblo0
      integer ibhi0
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0
      integer iweighthi0
      REAL*8 weight(
     &           iweightlo0:iweighthi0,
     &           0:nweightcomp-1)
      integer var
      integer icrse0
      integer ifine0
      integer ii0
      REAL*8 csum
      do var = 0, ncoarsecomp - 1
      do icrse0 = iblo0,ibhi0
            csum = 0
      do ii0 = iweightlo0,iweighthi0
               ifine0 = icrse0*ref_ratio + ii0 
               csum = csum + weight( ii0, 0) *
     &                 fine( ifine0, var )
      enddo
            coarse( icrse0, var ) = csum
      enddo
      end do
      return
      end
      subroutine NODEAVERAGEPOINT(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,pcrse
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0
     &           ,iweighthi0
     &           ,nweightcomp
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
      integer pcrse(0:0)
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0
      integer iweighthi0
      REAL*8 weight(
     &           iweightlo0:iweighthi0,
     &           0:nweightcomp-1)
      integer var
      integer ifine0
      integer ii0
      REAL*8 csum, weightpt, finept
      do var = 0, ncoarsecomp - 1
         csum = 0
      do ii0 = iweightlo0,iweighthi0
            ifine0 = pcrse(0)*ref_ratio + ii0 
            weightpt = weight( ii0, 0)
            finept = fine( ifine0, var )
            csum = csum + weightpt*finept
      enddo
         coarse(pcrse(0), var ) = csum
      end do
      return
      end
      subroutine NODEAVERAGE_GETWEIGHTS(
     &           weight
     &           ,iweightlo0
     &           ,iweighthi0
     &           ,nweightcomp
     &           ,ref_ratio
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nweightcomp
      integer iweightlo0
      integer iweighthi0
      REAL*8 weight(
     &           iweightlo0:iweighthi0,
     &           0:nweightcomp-1)
      integer ref_ratio
      integer ext, nxtrm
      integer ii0
      REAL*8 ref_scale
      ext = ref_ratio / 2
      ref_scale = (1.0d0) / (ref_ratio**1)
      do ii0 = iweightlo0,iweighthi0
         nxtrm = 0
         if (iabs(ii0) .eq. ext) nxtrm = nxtrm + 1 
         weight( ii0, 0) = ref_scale * (0.500d0)**nxtrm
      enddo
      return
      end
