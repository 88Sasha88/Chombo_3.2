#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine EXTRAPTOGHOSTVEL(
     &           vel
     &           ,ivello0
     &           ,ivelhi0
     &           ,nvelcomp
     &           ,iloboxlo0
     &           ,iloboxhi0
     &           ,haslo
     &           ,ihiboxlo0
     &           ,ihiboxhi0
     &           ,hashi
     &           ,facedir
     &           ,ncomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nvelcomp
      integer ivello0
      integer ivelhi0
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           0:nvelcomp-1)
      integer iloboxlo0
      integer iloboxhi0
      integer haslo
      integer ihiboxlo0
      integer ihiboxhi0
      integer hashi
      integer facedir
      integer ncomp
      integer ii,i, icomp
      
      ii = chf_id(facedir, 0)
      if(hashi.eq. 1) then
         do icomp = 0, ncomp-1
            
      do i = ihiboxlo0,ihiboxhi0

            vel(i,icomp) =
     $           2*vel(i-  ii, icomp)
     $           - vel(i-2*ii, icomp) 
            
      enddo
         enddo
      endif
      if(haslo .eq. 1) then
         do icomp = 0, ncomp-1
            
      do i = iloboxlo0,iloboxhi0

            vel(i,icomp) =
     $           2*vel(i+  ii, icomp)
     $           - vel(i+2*ii, icomp) 
            
      enddo
         enddo
      endif
      return 
      end
      subroutine APPLYOPVTOPNOBCS(
     &           lphfab
     &           ,ilphfablo0
     &           ,ilphfabhi0
     &           ,nlphfabcomp
     &           ,phifab
     &           ,iphifablo0
     &           ,iphifabhi0
     &           ,nphifabcomp
     &           ,acofab
     &           ,iacofablo0
     &           ,iacofabhi0
     &           ,eta0fab
     &           ,ieta0fablo0
     &           ,ieta0fabhi0
     &           ,eta1fab
     &           ,ieta1fablo0
     &           ,ieta1fabhi0
     &           ,eta2fab
     &           ,ieta2fablo0
     &           ,ieta2fabhi0
     &           ,lam0fab
     &           ,ilam0fablo0
     &           ,ilam0fabhi0
     &           ,lam1fab
     &           ,ilam1fablo0
     &           ,ilam1fabhi0
     &           ,lam2fab
     &           ,ilam2fablo0
     &           ,ilam2fabhi0
     &           ,dx
     &           ,alpha
     &           ,beta
     &           ,iregionlo0
     &           ,iregionhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlphfabcomp
      integer ilphfablo0
      integer ilphfabhi0
      REAL_T lphfab(
     &           ilphfablo0:ilphfabhi0,
     &           0:nlphfabcomp-1)
      integer nphifabcomp
      integer iphifablo0
      integer iphifabhi0
      REAL_T phifab(
     &           iphifablo0:iphifabhi0,
     &           0:nphifabcomp-1)
      integer iacofablo0
      integer iacofabhi0
      REAL_T acofab(
     &           iacofablo0:iacofabhi0)
      integer ieta0fablo0
      integer ieta0fabhi0
      REAL_T eta0fab(
     &           ieta0fablo0:ieta0fabhi0)
      integer ieta1fablo0
      integer ieta1fabhi0
      REAL_T eta1fab(
     &           ieta1fablo0:ieta1fabhi0)
      integer ieta2fablo0
      integer ieta2fabhi0
      REAL_T eta2fab(
     &           ieta2fablo0:ieta2fabhi0)
      integer ilam0fablo0
      integer ilam0fabhi0
      REAL_T lam0fab(
     &           ilam0fablo0:ilam0fabhi0)
      integer ilam1fablo0
      integer ilam1fabhi0
      REAL_T lam1fab(
     &           ilam1fablo0:ilam1fabhi0)
      integer ilam2fablo0
      integer ilam2fabhi0
      REAL_T lam2fab(
     &           ilam2fablo0:ilam2fabhi0)
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      integer iregionlo0
      integer iregionhi0
      real_t divf(0:CH_SPACEDIM-1)
      real_t aphi(0:CH_SPACEDIM-1)
      real_t  etaL(0:CH_SPACEDIM-1)
      real_t  etaH(0:CH_SPACEDIM-1)
      real_t  lamL(0:CH_SPACEDIM-1)
      real_t  lamH(0:CH_SPACEDIM-1)
      real_t  fluxL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  fluxH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiL(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  gphiH(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divuL(0:CH_SPACEDIM-1)
      real_t  divuH(0:CH_SPACEDIM-1)
      integer i, facedir, derivdir, veldir
      integer iif
      integer iid
      
      do i = iregionlo0,iregionhi0

      
      etaL(0) =eta0fab(i   )
                                  
      etaH(0) =eta0fab(i+1 )
      
      lamL(0) =lam0fab(i   )
                                  
      lamH(0) =lam0fab(i+1 )
      do facedir = 0, CH_SPACEDIM-1
         
         iif = chf_id(facedir, 0)
         do derivdir = 0, CH_SPACEDIM-1
            
            iid = chf_id(derivdir, 0)
            do veldir = 0, CH_SPACEDIM-1
               if(facedir .eq. derivdir) then
                  gphiH(veldir, derivdir ,facedir) = (phifab(i+iid,veldir) - phifab(i    ,veldir))/dx
                  gphiL(veldir, derivdir ,facedir) = (phifab(i    ,veldir) - phifab(i-iid,veldir))/dx
               else
                  gphiH(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phifab(i+iid+iif,veldir) - phifab(i-iid+iif,veldir) + 
     $                 phifab(i+iid    ,veldir) - phifab(i-iid    ,veldir)  )
                  gphiL(veldir, derivdir ,facedir) = (one/(four*dx))*(
     $                 phifab(i+iid-iif,veldir) - phifab(i-iid-iif,veldir) + 
     $                 phifab(i+iid    ,veldir) - phifab(i-iid    ,veldir)  )
               endif
            enddo
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         divuL(facedir) = zero
         divuH(facedir) = zero
         do veldir = 0, CH_SPACEDIM-1
            divuL(facedir)= divuL(facedir) + gphiL(veldir, veldir, facedir)
            divuH(facedir)= divuH(facedir) + gphiH(veldir, veldir, facedir)
         enddo
      enddo
      do facedir = 0, CH_SPACEDIM-1
         do veldir = 0, CH_SPACEDIM-1
            fluxL(veldir, facedir) = etaL(facedir)*(gphiL(facedir, veldir, facedir) + gphiL(veldir, facedir, facedir))
            fluxH(veldir, facedir) = etaH(facedir)*(gphiH(facedir, veldir, facedir) + gphiH(veldir, facedir, facedir))
            if(veldir .eq. facedir) then
               fluxL(veldir, facedir) = fluxL(veldir, facedir) + lamL(facedir)*divuL(facedir)
               fluxH(veldir, facedir) = fluxH(veldir, facedir) + lamH(facedir)*divuH(facedir)
            endif
         enddo
      enddo
      do veldir = 0, CH_SPACEDIM-1
         divf(veldir) = zero
         do facedir = 0, CH_SPACEDIM-1
            divf(veldir) = divf(veldir)
     $           + (fluxH(veldir, facedir)-fluxL(veldir, facedir))/dx
         enddo
      enddo
      
      aphi(0) = acofab(i)*phifab(i, 0)
      
      lphfab(i,0) =  alpha*aphi(0) + beta*divf(0)
      
      enddo
      ch_flops=ch_flops+(iregionhi0- iregionlo0+1)*(4 + 6*CH_SPACEDIM + 12*(CH_SPACEDIM*CH_SPACEDIM))
      return
      end
      subroutine GETFLUXVTOPNOBCS(
     &           flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,phifab
     &           ,iphifablo0
     &           ,iphifabhi0
     &           ,nphifabcomp
     &           ,etafab
     &           ,ietafablo0
     &           ,ietafabhi0
     &           ,lamfab
     &           ,ilamfablo0
     &           ,ilamfabhi0
     &           ,dx
     &           ,facedir
     &           ,beta
     &           ,ifaceboxlo0
     &           ,ifaceboxhi0
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           0:nfluxcomp-1)
      integer nphifabcomp
      integer iphifablo0
      integer iphifabhi0
      REAL_T phifab(
     &           iphifablo0:iphifabhi0,
     &           0:nphifabcomp-1)
      integer ietafablo0
      integer ietafabhi0
      REAL_T etafab(
     &           ietafablo0:ietafabhi0)
      integer ilamfablo0
      integer ilamfabhi0
      REAL_T lamfab(
     &           ilamfablo0:ilamfabhi0)
      REAL_T dx
      integer facedir
      REAL_T beta
      integer ifaceboxlo0
      integer ifaceboxhi0
      real_t  gphi(0:CH_SPACEDIM-1,0:CH_SPACEDIM-1)
      real_t  divu, eta, lam, fluxpt
      integer i , derivdir, veldir
      integer iif
      integer iid
      
      do i = ifaceboxlo0,ifaceboxhi0

      eta = etafab(i   );
      lam = lamfab(i   );
      
      iif = chf_id(facedir, 0)
      do derivdir = 0, CH_SPACEDIM-1
         
         iid = chf_id(derivdir, 0)
         do veldir = 0, CH_SPACEDIM-1
            if(facedir .eq. derivdir) then
               gphi(veldir, derivdir) = (phifab(i,veldir) - phifab(i-iid,veldir))/dx
            else
               gphi(veldir, derivdir) = (one/(four*dx))*(
     $              phifab(i+iid-iif,veldir) - phifab(i-iid-iif,veldir) + 
     $              phifab(i+iid    ,veldir) - phifab(i-iid    ,veldir)  )
            endif
         enddo
      enddo
      divu = zero
      do veldir = 0, CH_SPACEDIM-1
         divu = divu + gphi(veldir, veldir)
      enddo
      do veldir = 0, CH_SPACEDIM-1
         fluxpt = eta*(gphi(facedir, veldir) + gphi(veldir, facedir))
         if(veldir .eq. facedir) then
            fluxpt = fluxpt  + lam*divu
         endif
         flux(i, veldir) = beta*fluxpt
      enddo
      
      enddo
      return
      end
      subroutine PROLONGVTOP(
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
      integer ncomp, n
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
      subroutine RESTRICTRESVTOP(
     &           res
     &           ,ireslo0
     &           ,ireshi0
     &           ,nrescomp
     &           ,resfine
     &           ,iresfinelo0
     &           ,iresfinehi0
     &           ,nresfinecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,ncomp
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
      integer nresfinecomp
      integer iresfinelo0
      integer iresfinehi0
      REAL_T resfine(
     &           iresfinelo0:iresfinehi0,
     &           0:nresfinecomp-1)
      integer iregionlo0
      integer iregionhi0
      integer ncomp
      real_t denom
      integer n
      integer i
      integer ii
      denom = two
      do n = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

         
         ii = (i-abs(mod(i,2)))/2
         res(ii,n) = res(ii,n) + resfine(i,n)/denom
         
      enddo
      enddo
      return
      end
      subroutine CELLGRADVTOP(
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
      subroutine ADDGRADTOFLUXVTOP(
     &           flux
     &           ,ifluxlo0
     &           ,ifluxhi0
     &           ,nfluxcomp
     &           ,eta
     &           ,ietalo0
     &           ,ietahi0
     &           ,fluxcomp
     &           ,grad
     &           ,igradlo0
     &           ,igradhi0
     &           ,ngradcomp
     &           ,gradcomp
     &           ,gradtran
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
      integer ietalo0
      integer ietahi0
      REAL_T eta(
     &           ietalo0:ietahi0)
      integer fluxcomp
      integer ngradcomp
      integer igradlo0
      integer igradhi0
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           0:ngradcomp-1)
      integer gradcomp
      integer gradtran
      integer iregionfacelo0
      integer iregionfacehi0
      integer i
      real_t gradcontrib, trancontrib, etafac
      
      do i = iregionfacelo0,iregionfacehi0

      etafac      = eta(i)
      gradcontrib = grad(i, gradcomp)
      trancontrib = grad(i, gradtran)
      flux(i, fluxcomp) =  flux(i, fluxcomp) +
     $     eta(i)*
     $     ( grad(i, gradcomp)
     $     + grad(i, gradtran) )
      
      enddo
      return
      end
      subroutine GETFACEGRADVTOP(
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
      subroutine CELLDIVINCRVTOP(
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
      subroutine FACEDIVINCRVTOP(
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
      subroutine DECRINVRELCOEFVTOP(
     &           relcoef
     &           ,irelcoeflo0
     &           ,irelcoefhi0
     &           ,nrelcoefcomp
     &           ,eta
     &           ,ietalo0
     &           ,ietahi0
     &           ,netacomp
     &           ,lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,beta
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,dx
     &           ,idir
     &           ,ncomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nrelcoefcomp
      integer irelcoeflo0
      integer irelcoefhi0
      REAL_T relcoef(
     &           irelcoeflo0:irelcoefhi0,
     &           0:nrelcoefcomp-1)
      integer netacomp
      integer ietalo0
      integer ietahi0
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           0:netacomp-1)
      integer nlambdacomp
      integer ilambdalo0
      integer ilambdahi0
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           0:nlambdacomp-1)
      REAL_T beta
      integer iboxlo0
      integer iboxhi0
      REAL_T dx
      integer idir
      integer ncomp
      integer ii
      integer i
      integer icomp
      real_t lamh, laml,etah, etal, relcoold, relconew, incr
      
      ii = chf_id(idir, 0)
      do icomp = 0, ncomp-1
         
      do i = iboxlo0,iboxhi0

         relcoold = relcoef(i,icomp)
         etah =    eta(i+ii,0)
         etal =    eta(i   ,0)
         lamh = lambda(i+ii,0)
         laml = lambda(i   ,0)
         if(icomp .eq. idir) then
            incr =
     $         beta*(
     $        two*eta(i+ii,0) +
     $        lambda (i+ii,0) +
     $        two*eta(i   ,0) +
     $        lambda (i   ,0)
     $        )/(dx*dx)
         else
            incr =
     $         beta*(
     $        eta(i+ii,0) +
     $        eta(i   ,0)
     $        )/(dx*dx)
         endif
         relcoef(i, icomp) = relcoef(i,icomp) - incr
         relconew = relcoef(i,icomp)
         
      enddo
      enddo
      return
      end
      subroutine INITIALIZERELAXCOEF(
     &           relcoef
     &           ,irelcoeflo0
     &           ,irelcoefhi0
     &           ,nrelcoefcomp
     &           ,acoef
     &           ,iacoeflo0
     &           ,iacoefhi0
     &           ,alpha
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,ncomp
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer nrelcoefcomp
      integer irelcoeflo0
      integer irelcoefhi0
      REAL_T relcoef(
     &           irelcoeflo0:irelcoefhi0,
     &           0:nrelcoefcomp-1)
      integer iacoeflo0
      integer iacoefhi0
      REAL_T acoef(
     &           iacoeflo0:iacoefhi0)
      REAL_T alpha
      integer iboxlo0
      integer iboxhi0
      integer ncomp
      integer i
      integer icomp
      do icomp = 0, ncomp-1
         
      do i = iboxlo0,iboxhi0

         relcoef(i, icomp) = alpha*acoef(i)
         
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAVTOP(
     &           lambda
     &           ,ilambdalo0
     &           ,ilambdahi0
     &           ,nlambdacomp
     &           ,safety
     &           ,iboxlo0
     &           ,iboxhi0
     &           ,ncomp
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
      integer ncomp
      integer i
      integer icomp,jcomp
      real_t zeroval
      zeroval = 1.0e-20
         
      do i = iboxlo0,iboxhi0

         do icomp = 0, ncomp-1
            lambda(i, icomp) = safety/lambda(i, icomp)
         enddo
         
      enddo
      return
      end
      subroutine GSRBVTOP(
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
     &           ,ncomp
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
      integer ncomp
      integer i
      integer icomp
      REAL_T phio, lamo, rhso, lphio
      do icomp = 0, ncomp-1
         
      do i = icoloredboxlo0,icoloredboxhi0,2

         phio = phi(   i,icomp)
         lamo = lambda(i,icomp)
         rhso = rhs(   i,icomp)
         lphio = lphi(  i,icomp)
         phi(i, icomp) =
     $        phi(   i,icomp) +
     &        lambda(i,icomp)*(
     $        rhs(   i,icomp) -
     $        lphi(  i,icomp))
         
      enddo
      enddo
      return
      end
      subroutine ADDDIVFLUXDIRVTOP(
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
     &           ,ncomp
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
      integer ncomp
      integer facedir
      integer ii,i
      integer icomp
      
      ii = chf_id(facedir, 0)
      do icomp = 0, ncomp-1
         
      do i = iregionlo0,iregionhi0

         lhs(i, icomp) = lhs(i, icomp)
     $        +
     $        (flux(i+ii, icomp)
     $        -flux(i   , icomp))/dx
         
      enddo
      enddo
      return
      end
      subroutine LINEAREXTRAPVTOP(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,ighostBoxlo0
     &           ,ighostBoxhi0
     &           ,dir
     &           ,hiLo
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
      integer ighostBoxlo0
      integer ighostBoxhi0
      integer dir
      integer hiLo
      integer n
      integer i0
      integer ii0
      
      ii0= hiLo*CHF_ID(0, dir)

      do n=0, nphicomp-1
         
      do i0 = ighostBoxlo0,ighostBoxhi0

           phi(i0,n) = two*phi(i0-ii0,n)
     &                              - phi(i0-2*ii0,n)
         
      enddo
      enddo
      return
      end
      subroutine SLOPESVTOP(
     &           slopes
     &           ,islopeslo0
     &           ,islopeshi0
     &           ,nslopescomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,icBoxlo0
     &           ,icBoxhi0
     &           ,dir
     &           )

      implicit none
      integer*8 ch_flops
      COMMON/ch_timer/ ch_flops
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopescomp
      integer islopeslo0
      integer islopeshi0
      REAL_T slopes(
     &           islopeslo0:islopeshi0,
     &           0:nslopescomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer icBoxlo0
      integer icBoxhi0
      integer dir
      integer n, i0, ii0
      
      ii0= 1*CHF_ID(0, dir)

      do n=0, ncoarsecomp-1
         
      do i0 = icBoxlo0,icBoxhi0

         slopes(i0,n) = 0.5*(coarse(i0+ii0,n)
     &                                 -coarse(i0-ii0,n))
         
      enddo
      enddo
      return
      end
