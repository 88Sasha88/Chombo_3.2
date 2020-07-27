#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine PROLONGROP(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,m
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer iregionlo0
      integer iregionhi0
      integer m
      INTEGER ncomp, n
      integer i
      integer ii
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do i = iregionlo0,iregionhi0

          
          ii = (i-abs(mod(i,m)))/m
          phi(i,n) =  phi(i,n) +
     &        coarse(ii,n)
         
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESROP(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,resFine
     &           ,iresFinelo0
     &           ,iresFinehi0
     &           ,nresFinecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL_T res(
     &           ireslo0:ireshi0,
     &           0:nrescomp-1)
      integer nresFinecomp
      integer iresFinelo0
      integer iresFinehi0
      REAL_T resFine(
     &           iresFinelo0:iresFinehi0,
     &           0:nresFinecomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T denom
      integer n
      integer i
      integer ii
      integer ncomp
      ncomp = nrescomp
      denom = two
      do n = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

         
         ii = (i-abs(mod(i,2)))/2
         res(ii,n) = res(ii,n) + resFine(i,n)/denom
         
      enddo
      enddo
      return
      end
      subroutine CELLGRADROP(
     &           grad
     &           ,igradlo0
     &           ,igradhi0
     &           ,vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,igridlo0
     &           ,igridhi0
     &           ,dx
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradlo0
      integer igradhi0
      REAL_T grad(
     &           igradlo0:igradhi0)
      integer ivello0
      integer ivelhi0
      REAL_T vel(
     &           ivello0:ivelhi0)
      integer igridlo0
      integer igridhi0
      REAL_T dx
      integer divdir
      integer ii,i
      
      ii = chf_id(divdir, 0)
      
      do i = igridlo0,igridhi0

      grad(i) =
     $     (    vel(i+ii)
     $     -    vel(i-ii) )/(two*dx)
      
      enddo
      return
      end
      subroutine ADDGRADTOFLUXROP(
     &           flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,fluxComp
     &           ,grad
     &           ,igradlo0
     &           ,igradhi0
     &           ,ngradcomp
     &           ,gradComp
     &           ,sign
     &           ,iregionfacelo0
     &           ,iregionfacehi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer fluxComp
      integer ngradcomp
      integer igradlo0
      integer igradhi0
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           0:ngradcomp-1)
      integer gradComp
      REAL_T sign
      integer iregionfacelo0
      integer iregionfacehi0
      integer i
      
      do i = iregionfacelo0,iregionfacehi0

      flux(i, fluxcomp) =  flux(i, fluxcomp) +
     $     ( sign*grad(i, gradcomp))
      
      enddo
      return
      end
      subroutine GETFACEGRADROP(
     &           gradvelface
     &           ,igradvelfacelo0
     &           ,igradvelfacehi0
     &           ,gradvelcell
     &           ,igradvelcelllo0
     &           ,igradvelcellhi0
     &           ,velcomp
     &           ,ivelcomplo0
     &           ,ivelcomphi0
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,icenterboxlo0
     &           ,icenterboxhi0
     &           ,iloboxlo0
     &           ,iloboxhi0
     &           ,haslo
     &           ,ihiboxlo0
     &           ,ihiboxhi0
     &           ,hashi
     &           ,dx
     &           ,facedir
     &           ,divdir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0
      integer igradvelfacehi0
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0)
      integer igradvelcelllo0
      integer igradvelcellhi0
      REAL_T gradvelcell(
     &           igradvelcelllo0:igradvelcellhi0)
      integer ivelcomplo0
      integer ivelcomphi0
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0)
      integer iregionlo0
      integer iregionhi0
      integer icenterboxlo0
      integer icenterboxhi0
      integer iloboxlo0
      integer iloboxhi0
      integer haslo
      integer ihiboxlo0
      integer ihiboxhi0
      integer hashi
      REAL_T dx
      integer facedir
      integer divdir
      integer ii,i
      
      ii = chf_id(facedir, 0)
      if (facedir .eq. divdir) then
         
      do i = iregionlo0,iregionhi0

         gradvelface(i) =
     $        ( velcomp(i   )
     $        - velcomp(i-ii) )/dx
         
      enddo
      else
         
      do i = icenterboxlo0,icenterboxhi0

         gradvelface(i) =
     $        ( gradvelcell(i   )
     $        + gradvelcell(i-ii) )/two
         
      enddo
         if(haslo .eq. 1) then
            
      do i = iloboxlo0,iloboxhi0

            gradvelface(i) =
     $           (three*gradvelcell(i   )
     $           -      gradvelcell(i+ii))/two
            
      enddo
         endif
         if(hashi .eq. 1) then
            
      do i = ihiboxlo0,ihiboxhi0

            gradvelface(i) =
     $           (three*gradvelcell(i-  ii)
     $           -      gradvelcell(i-2*ii))/two
            
      enddo
         endif
      endif
      return
      end
      subroutine CELLDIVINCRROP(
     &           divvel
     &           ,idivvello0
     &           ,idivvelhi0
     &           ,vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,nvelcomp
     &           ,dx
     &           ,divdir
     &           ,iregionlo0
     &           ,iregionhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idivvello0
      integer idivvelhi0
      REAL_T divvel(
     &           idivvello0:idivvelhi0)
      integer nvelcomp
      integer ivello0
      integer ivelhi0
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           0:nvelcomp-1)
      REAL_T dx
      integer divdir
      integer iregionlo0
      integer iregionhi0
      integer ii,i
      
      ii = chf_id(divdir, 0)
      
      do i = iregionlo0,iregionhi0

      divvel(i) = divvel(i) +
     $     (    vel(i+ii,divdir)
     $     -    vel(i-ii,divdir) )/(two*dx)
      
      enddo
      return
      end
      subroutine FACEDIVINCRROP(
     &           divvel
     &           ,idivvello0
     &           ,idivvelhi0
     &           ,vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,nvelcomp
     &           ,gradvel
     &           ,igradvello0
     &           ,igradvelhi0
     &           ,ngradvelcomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,icenterboxlo0
     &           ,icenterboxhi0
     &           ,iloboxlo0
     &           ,iloboxhi0
     &           ,haslo
     &           ,ihiboxlo0
     &           ,ihiboxhi0
     &           ,hashi
     &           ,dx
     &           ,facedir
     &           ,divdir
     &           ,gradcomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idivvello0
      integer idivvelhi0
      REAL_T divvel(
     &           idivvello0:idivvelhi0)
      integer nvelcomp
      integer ivello0
      integer ivelhi0
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           0:nvelcomp-1)
      integer ngradvelcomp
      integer igradvello0
      integer igradvelhi0
      REAL_T gradvel(
     &           igradvello0:igradvelhi0,
     &           0:ngradvelcomp-1)
      integer iregionlo0
      integer iregionhi0
      integer icenterboxlo0
      integer icenterboxhi0
      integer iloboxlo0
      integer iloboxhi0
      integer haslo
      integer ihiboxlo0
      integer ihiboxhi0
      integer hashi
      REAL_T dx
      integer facedir
      integer divdir
      integer gradcomp
      integer ii,i
      
      ii = chf_id(facedir, 0)
      if (facedir .eq. divdir) then
         
      do i = iregionlo0,iregionhi0

         divvel(i) = divvel(i) +
     $        (    vel(i   ,facedir)
     $        -    vel(i-ii,facedir) )/dx
         
      enddo
      else
         
      do i = icenterboxlo0,icenterboxhi0

         divvel(i) = divvel(i) +
     $        ( gradvel(i   , gradcomp)
     $        + gradvel(i-ii, gradcomp) )/two
         
      enddo
         if(haslo .eq. 1) then
            
      do i = iloboxlo0,iloboxhi0

            divvel(i) = divvel(i) +
     $           (three*gradvel(i   , gradcomp)
     $           -      gradvel(i+ii, gradcomp))/two
            
      enddo
         endif
         if(hashi .eq. 1) then
            
      do i = ihiboxlo0,ihiboxhi0

            divvel(i) = divvel(i) +
     $           (three*gradvel(i-  ii, gradcomp)
     $           -      gradvel(i-2*ii, gradcomp))/two
            
      enddo
         endif
      endif
      return
      end
      subroutine DECRINVLAMBDAROP(
     &           lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,eta
     &           ,ietalo0
     &           ,ietahi0
     &           ,netacomp
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,beta
     &           ,dx
     &           ,idir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer netacomp
      integer ietalo0
      integer ietahi0
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           0:netacomp-1)
      integer iboxlo0
      integer iboxhi0
      REAL_T beta
      REAL_T dx
      integer idir
      integer ii
      integer i
      integer icomp
      REAL_T etahi, etalo
      integer ncomp
      ncomp = nlambdacomp
      
      ii = CHF_ID(idir, 0)
      do icomp = 0, ncomp-1
         
      do i = iboxlo0,iboxhi0

         etahi =  eta(i+ii,0)
         etalo =  eta(i   ,0)
         lambda(i, icomp) = lambda(i,icomp)
     $        - beta*(
     $        eta(i+ii,0) +
     $        eta(i   ,0))/(dx*dx)
         
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAROP(
     &           lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,safety
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      REAL_T safety
      integer iboxlo0
      integer iboxhi0
      integer i
      integer icomp
      integer ncomp
      ncomp = nlambdacomp
      do icomp = 0, ncomp-1
         
      do i = iboxlo0,iboxhi0

         lambda(i, icomp) =
     $        safety/(lambda(i,icomp))
         
      enddo
      enddo
      return
      end
      subroutine GSRBROP(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,lphi
     &           ,ilphilo0
     &           ,ilphihi0
     &           ,nlphicomp
     &           ,rhs
     &           ,irhslo0
     &           ,irhshi0
     &           ,nrhscomp
     &           ,lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,icoloredboxlo0
     &           ,icoloredboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer nlphicomp
      integer ilphilo0
      integer ilphihi0
      REAL_T lphi(
     &           ilphilo0:ilphihi0,
     &           0:nlphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           0:nrhscomp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      integer icoloredboxlo0
      integer icoloredboxhi0
      integer i
      integer icomp
      integer ncomp
      ncomp = nphicomp
      do icomp = 0, ncomp-1
         
      do i = icoloredBoxlo0,icoloredBoxhi0,2

         phi(i, icomp) =
     $        phi(   i,icomp) +
     &        lambda(i,icomp)*(
     $        rhs(   i,icomp) -
     $        lphi(  i,icomp))
         
      enddo
      enddo
      return
      end
      subroutine ADDDIVFLUXDIRROP(
     &           lhs
     &           ,ilhslo0
     &           ,ilhshi0
     &           ,nlhscomp
     &           ,flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,dx
     &           ,facedir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlhscomp
      integer ilhslo0
      integer ilhshi0
      REAL_T lhs(
     &           ilhslo0:ilhshi0,
     &           0:nlhscomp-1)
      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T dx
      integer facedir
      integer ii,i
      integer icomp
      integer ncomp
      ncomp = nlhscomp
      
      ii = CHF_ID(facedir, 0)
      do icomp = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

         lhs(i, icomp) = lhs(i, icomp) +
     $        (flux(i+ii, icomp)
     $        -flux(i   , icomp))/dx
         
      enddo
      enddo
      return
      end
