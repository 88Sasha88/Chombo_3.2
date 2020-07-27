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
      REAL*8 rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           0:nrhscomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer idomainlo0,idomainlo1
      integer idomainhi0,idomainhi1
      REAL*8 dx
      REAL*8 rhono
      REAL*8 rno
      integer iprob
      integer i,j
      REAL*8 xi,yi
      REAL*8 rad
      REAL*8 xloc, xlo,xhi
      integer nv
      REAL*8 rhsexact, signComp
      xlo = idomainlo0*dx
      xhi =(idomainhi0+1)*dx
      xloc = (0.500d0)*(xlo + xhi)
      signComp = -(1.0d0)
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
            xi = (2.0d0)*(3.14159265358979323846264338327950288d0)*dx*(i
     &+0.5)
            yi = (2.0d0)*(3.14159265358979323846264338327950288d0)*dx*(j
     &+0.5)
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
      REAL*8 rhono
      REAL*8 rno
      REAL*8 rad
      integer iprob
      REAL*8 rhsexact
      REAL*8  rhs, radrat
      radrat = rad/rno
      if(iprob.eq.0) then
         if(rad .lt. rno) then
            rhs = rhono
         else
            rhs = (0.0d0)
         endif
      elseif(iprob.eq.1) then
         if(rad .lt. rno) then
            rhs = rhono*((2.0d0)*radrat**3 - (3.0d0)*radrat**2 + 1)
         else
            rhs = (0.0d0)
         endif
      else
         call MAYDAY_ERROR()
      endif
      rhsexact = rhs
      return
      end
