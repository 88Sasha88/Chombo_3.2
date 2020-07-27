#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FACENODEBC(
     &           state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,neumfac
     &           ,ineumfaclo0
     &           ,ineumfachi0
     &           ,nneumfaccomp
     &           ,dircfac
     &           ,idircfaclo0
     &           ,idircfachi0
     &           ,ndircfaccomp
     &           ,inhmval
     &           ,iinhmvallo0
     &           ,iinhmvalhi0
     &           ,ninhmvalcomp
     &           ,ifaceboxlo0
     &           ,ifaceboxhi0
     &           ,idir
     &           ,side
     &           ,dx
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer nneumfaccomp
      integer ineumfaclo0
      integer ineumfachi0
      REAL_T neumfac(
     &           ineumfaclo0:ineumfachi0,
     &           0:nneumfaccomp-1)
      integer ndircfaccomp
      integer idircfaclo0
      integer idircfachi0
      REAL_T dircfac(
     &           idircfaclo0:idircfachi0,
     &           0:ndircfaccomp-1)
      integer ninhmvalcomp
      integer iinhmvallo0
      integer iinhmvalhi0
      REAL_T inhmval(
     &           iinhmvallo0:iinhmvalhi0,
     &           0:ninhmvalcomp-1)
      integer ifaceboxlo0
      integer ifaceboxhi0
      integer idir
      integer side
      REAL_T dx
      integer startcomp
      integer endcomp
      REAL_T nfac, dfac, ival, sval,denom,numer
      integer ncomp,nc
      integer i0, ii0
      ncomp = nstatecomp
      if(ncomp .ne. nneumfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ndircfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ninhmvalcomp) then
          call MAYDAY_ERROR()
      endif
      if ((side .ne. -1) .and. (side .ne. 1)) then
          call MAYDAY_ERROR()
      endif
      
      ii0= side*CHF_ID(0, idir)

      do nc = startcomp, endcomp
          
      do i0 = ifaceboxlo0,ifaceboxhi0

              nfac = neumfac(i0, nc)
              dfac = dircfac(i0, nc)
              ival = inhmval(i0, nc)
              sval = state(i0-ii0, nc)
              denom = dfac + side*(nfac/dx)
              numer = ival + side*(nfac/dx)*sval
#ifndef NDEBUG
              if (abs(denom) .lt. 1.0e-9) then
                  call MAYDAY_ERROR()
              endif
#endif
              state(i0, nc) = numer/denom
          
      enddo
      enddo
      return
      end
