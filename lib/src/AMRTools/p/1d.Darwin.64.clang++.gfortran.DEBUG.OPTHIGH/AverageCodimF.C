#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine AVERAGECODIM(
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
     &           ,refFactor
     &           ,refDim
     &           ,ibreflo0
     &           ,ibrefhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer iboxlo0
      integer iboxhi0
      integer refRatio
      integer refFactor
      integer refDim
      integer ibreflo0
      integer ibrefhi0
      integer var
      integer ic0
      integer ip0
      integer ii0
      real_t refScale,coarseSum
      refScale = one/(refFactor**refDim)
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iboxlo0,iboxhi0

            
               ip0 = ic0*refRatio
            coarseSum = 0
            
      do ii0 = ibreflo0,ibrefhi0

               coarseSum = coarseSum + fine(ip0+ii0,var)
            
      enddo
            coarse(ic0,var) = coarseSum*refScale
         
      enddo
      enddo
      return
      end
      subroutine AVERAGECODIMHARMONIC(
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
     &           ,refFactor
     &           ,refDim
     &           ,ibreflo0
     &           ,ibrefhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer iboxlo0
      integer iboxhi0
      integer refRatio
      integer refFactor
      integer refDim
      integer ibreflo0
      integer ibrefhi0
      integer var
      integer ic0
      integer ip0
      integer ii0
      real_t refScale,coarseSum
      refScale = one/(refFactor**refDim)
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iboxlo0,iboxhi0

            
               ip0 = ic0*refRatio
            coarseSum = 0
            
      do ii0 = ibreflo0,ibrefhi0

               coarseSum = coarseSum + one/
     &            fine(ip0+ii0,var)
            
      enddo
            coarse(ic0,var) = one/(coarseSum*refScale)
         
      enddo
      enddo
      return
      end
