      subroutine REGCORRECTTANVEL(
     &           vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,grad
     &           ,igradlo0
     &           ,igradhi0
     &           ,iinteriorboxlo0
     &           ,iinteriorboxhi0
     &           ,veldir
     &           ,graddir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivello0
      integer ivelhi0
      REAL*8 vel(
     &           ivello0:ivelhi0)
      integer igradlo0
      integer igradhi0
      REAL*8 grad(
     &           igradlo0:igradhi0)
      integer iinteriorboxlo0
      integer iinteriorboxhi0
      integer veldir
      integer graddir
      integer ioffvel 
      integer ioffgrad
      integer i
      REAL*8 correction, factor
      ioffvel = chf_id(0,veldir)
      ioffgrad = chf_id(0,graddir)
      factor = (1.0d0)/(4.0d0)
      do i = iinteriorboxlo0,iinteriorboxhi0
      correction = factor*
     $     ( grad(i                 )
     $     + grad(i+ioffgrad-ioffvel)
     $     + grad(i         -ioffvel)
     $     + grad(i+ioffgrad        ))
      vel(i) = vel(i) - correction
      enddo
      ch_flops=ch_flops+(iinteriorboxhi0- iinteriorboxlo0+1)*5+1
      return
      end
      subroutine REGCORRECTTANVELVARDENS(
     &           vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,grad
     &           ,igradlo0
     &           ,igradhi0
     &           ,oneoverdens
     &           ,ioneoverdenslo0
     &           ,ioneoverdenshi0
     &           ,iinteriorboxlo0
     &           ,iinteriorboxhi0
     &           ,veldir
     &           ,graddir
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivello0
      integer ivelhi0
      REAL*8 vel(
     &           ivello0:ivelhi0)
      integer igradlo0
      integer igradhi0
      REAL*8 grad(
     &           igradlo0:igradhi0)
      integer ioneoverdenslo0
      integer ioneoverdenshi0
      REAL*8 oneoverdens(
     &           ioneoverdenslo0:ioneoverdenshi0)
      integer iinteriorboxlo0
      integer iinteriorboxhi0
      integer veldir
      integer graddir
      integer ioffvel 
      integer ioffgrad
      integer i
      REAL*8 correction, factor
      ioffvel = chf_id(0,veldir)
      ioffgrad = chf_id(0,graddir)
      factor = (1.0d0)/(4.0d0)
      do i = iinteriorboxlo0,iinteriorboxhi0
      correction = factor*
     $     ( grad(i                 )
     $     + grad(i+ioffgrad-ioffvel)
     $     + grad(i         -ioffvel)
     $     + grad(i+ioffgrad        ))
      vel(i) = vel(i) - correction*oneoverdens(i)
      enddo
      ch_flops=ch_flops+(iinteriorboxhi0- iinteriorboxlo0+1)*5+1
      return
      end
      subroutine MACDIVERGEF(
     &           idcalclo0
     &           ,idcalchi0
     &           ,divf
     &           ,idivflo0
     &           ,idivfhi0
     &           ,ndivfcomp
     &           ,flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,facedir
     &           ,nconserved
     &           ,dx
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0
      integer idcalchi0
      integer ndivfcomp
      integer idivflo0
      integer idivfhi0
      REAL*8 divf(
     &           idivflo0:idivfhi0,
     &           0:ndivfcomp-1)
      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL*8 flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dx
      integer i
      integer ioff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      spacedim = 1
      do iv = 0,nconserved - 1
      do i = idcalclo0,idcalchi0
         divf(i,iv) = divf(i,iv) +
     &        (flux(i+ioff,iv)
     &        -flux(i     ,iv))/dx
      enddo
      enddo
      ch_flops=ch_flops+(idcalchi0- idcalclo0+1)*3
      return
      end
      subroutine MACGRADPHI(
     &           gradphi
     &           ,igradphilo0
     &           ,igradphihi0
     &           ,phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,facedir
     &           ,dx
     &           ,ifaceboxlo0
     &           ,ifaceboxhi0
     &           )
      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradphilo0
      integer igradphihi0
      REAL*8 gradphi(
     &           igradphilo0:igradphihi0)
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     &           iphilo0:iphihi0)
      integer facedir
      REAL*8 dx
      integer ifaceboxlo0
      integer ifaceboxhi0
      integer i
      integer ioff
      ioff = chf_id(0,facedir)
      do i = ifaceboxlo0,ifaceboxhi0
      gradphi(i) =
     &     ( phi(i     )
     &     - phi(i-ioff)
     &     )/dx
      enddo
      ch_flops=ch_flops+(ifaceboxhi0- ifaceboxlo0+1)*2
      return
      end
