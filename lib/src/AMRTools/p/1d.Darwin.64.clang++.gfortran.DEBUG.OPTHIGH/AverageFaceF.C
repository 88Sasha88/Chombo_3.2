#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine AVERAGEFACE(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,icrseBoxlo0
     &           ,icrseBoxhi0
     &           ,dir
     &           ,nRef
     &           ,refFactor
     &           ,irefBoxlo0
     &           ,irefBoxhi0
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
      integer icrseBoxlo0
      integer icrseBoxhi0
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0
      integer irefBoxhi0
      integer ic0
      integer ifine0
      integer var
      integer ii0
      REAL_T crseSum, ref_scale
      ref_scale = (one/refFactor)**(CH_SPACEDIM-1)
      do var=0, ncoarsecomp-1
         
      do ic0 = icrseBoxlo0,icrseBoxhi0

         crseSum = 0
         
      do ii0 = irefBoxlo0,irefBoxhi0

         
         ifine0=ic0*nRef+ii0
            crseSum = crseSum + fine(ifine0,var)
            
      enddo
            coarse(ic0,var) = ref_scale*crseSum
          
      enddo
       enddo
       return
       end
      subroutine AVERAGEFACEHARMONIC(
     &           coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,icrseBoxlo0
     &           ,icrseBoxhi0
     &           ,dir
     &           ,nRef
     &           ,refFactor
     &           ,irefBoxlo0
     &           ,irefBoxhi0
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
      integer icrseBoxlo0
      integer icrseBoxhi0
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0
      integer irefBoxhi0
      integer ic0
      integer ifine0
      integer var
      integer ii0
      REAL_T crseSum, ref_scale
      ref_scale = (one/refFactor)**(CH_SPACEDIM-1)
      do var=0, ncoarsecomp-1
         
      do ic0 = icrseBoxlo0,icrseBoxhi0

         crseSum = 0
         
      do ii0 = irefBoxlo0,irefBoxhi0

         
         ifine0=ic0*nRef+ii0
            crseSum = crseSum + one/fine(ifine0,var)
            
      enddo
            coarse(ic0,var) = one/(ref_scale*crseSum)
          
      enddo
       enddo
       return
       end
