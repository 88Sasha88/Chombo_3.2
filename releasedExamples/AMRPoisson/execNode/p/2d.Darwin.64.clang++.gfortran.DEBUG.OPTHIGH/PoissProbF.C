#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine GETRHSPOIS(
     &           rhs
     &           ,irhslo0,irhslo1
     &           ,irhshi0,irhshi1
     &           ,nrhscomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,idomainlo0,idomainlo1
     &           ,idomainhi0,idomainhi1
     &           ,dx
     &           ,rhono
     &           ,rno
     &           ,iprob
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrhscomp
      integer irhslo0,irhslo1
      integer irhshi0,irhshi1
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idomainlo0,idomainlo1
      integer idomainhi0,idomainhi1
      REAL_T dx
      REAL_T rhono
      REAL_T rno
      integer iprob
      integer i,j
      REAL_T xi,yi
      REAL_T rad
      REAL_T xloc, xlo,xhi
      integer nv
      REAL_T rhsexact, signComp
      xlo = idomainlo0*dx
      xhi =(idomainhi0+1)*dx
      xloc = half*(xlo + xhi)
      signComp = -one
      if (iprob.ne.2) then
         do nv = 0, nrhscomp - 1
            signComp = -signComp
            
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            
            xi = dx*(i+0.5) - xloc
            yi = dx*(j+0.5) - xloc
            rad = sqrt(
     $           xi*xi+ yi*yi
     $           )
            call getrhsexact(rhono, rno, rad, iprob, rhsexact)
            rhs(i,j,nv) =  signComp*rhsexact
            
      enddo
      enddo
         enddo
      else
         do nv = 0, nrhscomp - 1
            signComp = -signComp
            
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

            
            xi = two*Pi*dx*(i+0.5)
            yi = two*Pi*dx*(j+0.5)
            rhs(i,j,nv) =sin(xi)+ sin(yi)
            rhs(i,j,nv) = signComp*rhs(i,j,nv)
            
      enddo
      enddo
         enddo
      endif
      return
      end
      subroutine GETRHSEXACT(
     &           rhono
     &           ,rno
     &           ,rad
     &           ,iprob
     &           ,rhsexact
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      REAL_T rhono
      REAL_T rno
      REAL_T rad
      integer iprob
      REAL_T rhsexact
      REAL_T  rhs, radrat
      radrat = rad/rno
      if(iprob.eq.0) then
         if(rad .lt. rno) then
            rhs = rhono
         else
            rhs = zero
         endif
      elseif(iprob.eq.1) then
         if(rad .lt. rno) then
            rhs = rhono*(two*radrat**3 - three*radrat**2 + 1)
         else
            rhs = zero
         endif
      else
         call MAYDAY_ERROR()
      endif
      rhsexact = rhs
      return
      end
