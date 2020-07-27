#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine TRAPWEIGHTS(
     &           wt
     &           ,iwtlo0
     &           ,iwthi0
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dv
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer iwtlo0
      integer iwthi0
      REAL_T wt(
     &           iwtlo0:iwthi0)
      integer iregionlo0
      integer iregionhi0
      REAL_T dv
      integer i0
      
      do i0 = iregionlo0,iregionhi0

         wt(i0) = dv
         
         if ((i0 .eq. iregionlo0) .or.
     &        (i0 .eq. iregionhi0)) then
             wt(i0) = wt(i0) * half
         endif 
      
      enddo
      return
      end
      subroutine NODEDOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0
     &           ,iafabhi0
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0
     &           ,ibfabhi0
     &           ,nbfabcomp
     &           ,wfab
     &           ,iwfablo0
     &           ,iwfabhi0
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0
      integer iafabhi0
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0
      integer ibfabhi0
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           0:nbfabcomp-1)
      integer iwfablo0
      integer iwfabhi0
      REAL_T wfab(
     &           iwfablo0:iwfabhi0)
      integer iregionlo0
      integer iregionhi0
      integer startcomp
      integer endcomp
      integer i0
      integer nv
      dotprodout = 0
      do nv=startcomp,endcomp,1
         
      do i0 = iregionlo0,iregionhi0

            dotprodout = dotprodout +
     &        afab(i0,nv)*
     &        bfab(i0,nv)*
     &        wfab(i0)
         
      enddo
      enddo
      return
      end
      subroutine NODENORM(
     &           normout
     &           ,p
     &           ,afab
     &           ,iafablo0
     &           ,iafabhi0
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0
     &           ,iwfabhi0
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T normout
      integer p
      integer nafabcomp
      integer iafablo0
      integer iafabhi0
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           0:nafabcomp-1)
      integer iwfablo0
      integer iwfabhi0
      REAL_T wfab(
     &           iwfablo0:iwfabhi0)
      integer iregionlo0
      integer iregionhi0
      integer startcomp
      integer endcomp
      integer i0
      integer nv
      REAL_T pwrinv
      normout = 0
      if (p .eq. 1) then
         do nv = startcomp, endcomp, 1
            
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0) *
     &           abs(afab(i0, nv))
            
      enddo
         enddo
      elseif (p .eq. 2) then
         do nv = startcomp, endcomp, 1
            
      do i0 = iregionlo0,iregionhi0

            normout = normout +
     &           wfab(i0) *
     &           afab(i0, nv) *
     &           afab(i0, nv)
            
      enddo
         enddo
         normout = sqrt(normout)
      else
         do nv = startcomp, endcomp, 1
            
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0) *
     &           (afab(i0, nv)**p)
            
      enddo
         enddo
         pwrinv = one / p
         normout = normout**pwrinv
      endif
      return
      end
      subroutine NODEINTEGRAL(
     &           ans
     &           ,afab
     &           ,iafablo0
     &           ,iafabhi0
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0
     &           ,iwfabhi0
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T ans
      integer nafabcomp
      integer iafablo0
      integer iafabhi0
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           0:nafabcomp-1)
      integer iwfablo0
      integer iwfabhi0
      REAL_T wfab(
     &           iwfablo0:iwfabhi0)
      integer iregionlo0
      integer iregionhi0
      integer startcomp
      integer endcomp
      integer i0
      integer nv
      ans = 0
      do nv = startcomp, endcomp, 1
         
      do i0 = iregionlo0,iregionhi0

            ans = ans +
     &        wfab(i0) * afab(i0, nv)
         
      enddo
      enddo
      return
      end
      subroutine NODEMAXNORM(
     &           normout
     &           ,afab
     &           ,iafablo0
     &           ,iafabhi0
     &           ,nafabcomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T normout
      integer nafabcomp
      integer iafablo0
      integer iafabhi0
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           0:nafabcomp-1)
      integer iregionlo0
      integer iregionhi0
      integer startcomp
      integer endcomp
      integer i0
      integer nv
      REAL_T this
      normout = 0
      do nv = startcomp, endcomp, 1
         
      do i0 = iregionlo0,iregionhi0

            this = abs(afab(i0, nv))
            if (this .gt. normout) then
               normout = this
            endif
         
      enddo
      enddo
      return
      end
