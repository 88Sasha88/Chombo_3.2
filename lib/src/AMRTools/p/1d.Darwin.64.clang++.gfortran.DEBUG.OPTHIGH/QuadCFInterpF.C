#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine QUADINTERP(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,phistar
     &           ,iphistarlo0
     &           ,iphistarhi0
     &           ,nphistarcomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,ihilo
     &           ,h
     &           ,idir
     &           ,scomp
     &           ,ecomp
     &           ,nref
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nphistarcomp
      integer iphistarlo0
      integer iphistarhi0
      REAL_T phistar(
     &           iphistarlo0:iphistarhi0,
     &           0:nphistarcomp-1)
      integer iboxlo0
      integer iboxhi0
      integer ihilo
      REAL_T h
      integer idir
      integer scomp
      integer ecomp
      integer nref
      integer i0
      integer ii0
      integer n
      real_t x, pa, pb, ps, a, b, frac, denom, xsquared
      real_t  mult, invh
      frac = two/(h*h)
      denom =  nref*nref + 4*nref + 3
      mult = frac/denom
      invh = one / h
      x = two * h
      xsquared = four * h*h
      
      ii0= ihilo*CHF_ID(0, idir)

      do n=scomp,ecomp
         
      do i0 = iboxlo0,iboxhi0

            pa = phi(i0 -2*ii0,n)
            pb = phi(i0 -ii0,n)
            ps = phistar(i0 +ii0,n)
            a  = mult*(two*ps + (nref+1)*pa - (nref+3)*pb)
            b  = (pb-pa)*invh - a*h
            phi(i0,n) = xsquared*a + b*x + pa
         
      enddo
      enddo
      return
      end
      subroutine PHISTAR(
     &           fPhiStar
     &           ,ifPhiStarlo0
     &           ,ifPhiStarhi0
     &           ,nfPhiStarcomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,phic
     &           ,iphiclo0
     &           ,iphichi0
     &           ,nphiccomp
     &           ,coarslope
     &           ,icoarslopelo0
     &           ,icoarslopehi0
     &           ,ncoarslopecomp
     &           ,coarcurva
     &           ,icoarcurvalo0
     &           ,icoarcurvahi0
     &           ,ncoarcurvacomp
     &           ,coarmixed
     &           ,icoarmixedlo0
     &           ,icoarmixedhi0
     &           ,ncoarmixedcomp
     &           ,dxf
     &           ,ivar
     &           ,dir
     &           ,sign
     &           ,nRef
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfPhiStarcomp
      integer ifPhiStarlo0
      integer ifPhiStarhi0
      REAL_T fPhiStar(
     &           ifPhiStarlo0:ifPhiStarhi0,
     &           0:nfPhiStarcomp-1)
      integer iregionlo0
      integer iregionhi0
      integer nphiccomp
      integer iphiclo0
      integer iphichi0
      REAL_T phic(
     &           iphiclo0:iphichi0,
     &           0:nphiccomp-1)
      integer ncoarslopecomp
      integer icoarslopelo0
      integer icoarslopehi0
      REAL_T coarslope(
     &           icoarslopelo0:icoarslopehi0,
     &           0:ncoarslopecomp-1)
      integer ncoarcurvacomp
      integer icoarcurvalo0
      integer icoarcurvahi0
      REAL_T coarcurva(
     &           icoarcurvalo0:icoarcurvahi0,
     &           0:ncoarcurvacomp-1)
      integer ncoarmixedcomp
      integer icoarmixedlo0
      integer icoarmixedhi0
      REAL_T coarmixed(
     &           icoarmixedlo0:icoarmixedhi0,
     &           0:ncoarmixedcomp-1)
      REAL_T dxf
      integer ivar
      integer dir
      integer sign
      integer nRef
      REAL_T xf1, xc1, xf2, xc2, x1, x2, dxc
      REAL_T aa, update1, update2, update3
      integer i0
      integer ii0
      integer ir0
      integer ic(0:CH_SPACEDIM-1)
      integer ivf(0:CH_SPACEDIM-1)
      integer YOU(1:2, 0:2), you1, you2
      data YOU / 1, 2, 0, 2, 0, 1 /
#if CH_SPACEDIM > 3
      call MAYDAY_ERROR()
#else
      dxc = nRef * dxf
      you1 = YOU(1,dir)
      you2 = YOU(2,dir)
      
      ii0= sign*CHF_ID(0, dir)

      
      do ir0 = iregionlo0,iregionhi0

         
         ic(0)=ir0/nRef
         
         ivf(0)=ir0
         
         i0=ir0+ii0
         xf1 = (ivf(you1)+half)*dxf
         xc1 = (  ic(you1)+half)*dxc
         xf2 = (ivf(you2)+half)*dxf
         xc2 = (  ic(you2)+half)*dxc
         x1 = xf1-xc1
         x2 = xf2-xc2
         aa= phic(ic(0),ivar)
         update1=x1*coarslope(ic(0),you1) +
     &        half*x1*x1*coarcurva(ic(0),you1)
         update2=x2*coarslope(ic(0),you2) +
     &        half*x2*x2*coarcurva(ic(0),you2)
         update3=x1*x2*coarmixed(ic(0),0)
         fPhiStar(i0,ivar) = aa+update1+update2+update3
      
      enddo
#endif
      return
      end
