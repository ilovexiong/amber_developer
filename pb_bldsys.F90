! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map atomic charges onto grid
subroutine pb_crggrd ( natom,atmfirst,atmlast,xm,ym,zm,icrd,gcrg,chgrd )

   use poisson_boltzmann, only : ligand,multiblock,realflag,level,nfocus
   implicit none

   integer natom,atmfirst, atmlast
   integer xm, ym, zm
   integer icrd(3,natom)
   _REAL_ gcrg(8,natom)
   _REAL_ chgrd(xm,ym,zm)

   integer iatm, i, j, k
   do iatm = atmfirst, atmlast
      if ( ligand .or. multiblock ) then
         if ( level == nfocus .and. realflag(iatm) == 0 ) cycle
      endif
      !if ( level == nfocus .and. ligand     .and. realflag(iatm) == 0 ) cycle
      !if ( level == nfocus .and. multiblock .and. realflag(iatm) == 0 ) cycle
      i = icrd(1, iatm); j = icrd(2, iatm); k = icrd(3, iatm)
      if ( i<1 .or. j<1 .or. k<1 .or. i+1>xm .or. j+1>ym .or. k+1>zm ) then
         write(6,'(6i6)') i,j,k,xm,ym,zm
         write(6,'(a)') "pb_crggrd: index exceeds box range, bail."
         call mexit(6,1)
      endif
      chgrd(i  ,j  ,k  ) = chgrd(i  ,j  ,k  ) + gcrg(1, iatm)
      chgrd(i+1,j  ,k  ) = chgrd(i+1,j  ,k  ) + gcrg(2, iatm)
      chgrd(i  ,j+1,k  ) = chgrd(i  ,j+1,k  ) + gcrg(3, iatm)
      chgrd(i+1,j+1,k  ) = chgrd(i+1,j+1,k  ) + gcrg(4, iatm)
      chgrd(i  ,j  ,k+1) = chgrd(i  ,j  ,k+1) + gcrg(5, iatm)
      chgrd(i+1,j  ,k+1) = chgrd(i+1,j  ,k+1) + gcrg(6, iatm)
      chgrd(i  ,j+1,k+1) = chgrd(i  ,j+1,k+1) + gcrg(7, iatm)
      chgrd(i+1,j+1,k+1) = chgrd(i+1,j+1,k+1) + gcrg(8, iatm)
   end do

end subroutine pb_crggrd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular ion map assignment
subroutine pb_ionmap( pbverbose,ifcap,natom,memopt,mzmin,mzmax,iprob,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              outflag,gcrd,radi,atmsas,insas,saltgrd)

   implicit none

   ! passed variables
 
   logical pbverbose
   integer ifcap
   integer natom
   integer memopt
   _REAL_ mzmin, mzmax
   _REAL_ iprob, h, gox, goy, goz
   integer xm, ym, zm, xmymzm
   integer outflag(*)
   _REAL_ gcrd(3,*), radi(*)
   integer atmsas(*), insas(*)
   _REAL_ saltgrd(xmymzm)

   ! local variables

   integer iatm, xmymzm_ext
   _REAL_ rh, range0, range1, xi, yi, zi

   rh = 1.0d0/h

   ! local array setup

   !xmymzm_ext = xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8
   xmymzm_ext = (xm+2)*(ym+2)*(zm+2)

   ! resetting atmsas and insas, -4 refers to bulk water and salt.

   insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0

   ! mark grid points within the Stern layer as -3

   ! part a: mark the membrane portion, i.e. with mzmin < z < mzmax with iprob
   ! added

   if ( memopt > 0 ) then
      range0 = mzmin-iprob*rh
      range1 = mzmax+iprob*rh
      call exstslab( -3, insas )
   end if

   ! part b: mark each solute atom one by one

   do iatm = 1, natom
      if (ifcap == 5 .and. outflag(iatm) == 1) cycle
      range0 = radi(iatm)
      if ( range0 == 0.0d0 ) cycle
      range1 = (range0+iprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      call exstsph( iatm, -3, insas, atmsas )
   end do

   ! part c: set up the ion exclusion map

   saltgrd(1:xmymzm) = 1.0d0
   call ionmap( insas, saltgrd )

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Exlusion of ions from protein interior.
subroutine ionmap ( insas,saltgrd )

   ! Passed variables

   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ saltgrd(xm,ym,zm)

   ! Local variables

   integer i, j, k
   ! for InsightII display
   !_REAL_ g(3)

   ! for InsightII display
   !open (unit=55, file='ions.dot')
   !write (55, '("DOTS")')
   do k = 1, zm; do j = 1, ym; do i = 1, xm
      if ( insas(i,j,k) /= -4 ) then
         saltgrd(i,j,k) = 0.0d0
         ! for InsightII display
         !g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
         !write (55,'(4(f8.3,2x))') g(1:3), 300.
      end if
   end do; end do; end do
   ! for InsightII display
   !close(55)

end subroutine ionmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within an rectangular slab in parallel to the xy plane
subroutine exstslab( dielslab,inslab )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a slab (dielectric constant dielsph) between
   ! zmin=range0 and zmax=range1 (in the FD grid unit and frame) as dielpsh.
   ! Note the difference with the routine for atoms, which also saves atom ID for 
   ! potential force calculations.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer dielslab
   integer inslab(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer lowk, highk

   lowk = max(0,ceiling(range0)); highk = min(zm+1,floor(range1))

   inslab(0:xm+1,0:ym+1,lowk:highk) = dielslab

end subroutine exstslab             
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic stern spheres
subroutine exstsph( iatm,dielsph,insph,inatm )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! Passed variables

   integer iatm
   integer dielsph
   integer insph(0:xm+1,0:ym+1,0:zm+1)
   integer inatm(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer i, j, k
   integer lowi, lowj, lowk
   integer highi, highj, highk
   _REAL_ range2, range3!, d, d2, r

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               if ( insph(i,j,k) == dielsph ) cycle
               insph(i,j,k) = dielsph; inatm(i,j,k) = iatm
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
 
end subroutine exstsph

end subroutine pb_ionmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment
!     call pb_exmol( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
!             smoothopt,dprob,epsin,epsout,epsmem,membraneopt,mzmin,mzmax,&
!             h,gox,goy,goz,xm,ym,zm,xmymzm,&
!             level,nfocus,&
!             narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
!             outflag,gcrd,acrd,radi,radip3,&
!             marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
!             atmsas,insas,lvlset,mlvlset,zv,epsx,epsy,epsz,&
!             iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)
subroutine pb_exmol( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
              smoothopt,dprob,epsin,epsout,epsmem,memopt,mzmin,mzmax,&
              h,gox,goy,goz,xm,ym,zm,xmymzm,&
              level,nfocus,&
              narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
              outflag,gcrd,acrd,radi,radip3,&
              marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
              atmsas,insas,lvlset,mlvlset,zv,epsx,epsy,epsz,&
              iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)

   use poisson_boltzmann, only: ligand, multiblock, mthick,&
                              mctrdz, outlvlset, outmlvlset, poretype,&
                              poreradi, cutfd
   use density_surface, only: density_lvlset, irreg_init, irreg_grd, bndatm
   use pbtimer_module
   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap, ipb, natom, smoothopt
   _REAL_ dprob, epsin, epsout, epsmem
   integer memopt
   _REAL_ mzmin, mzmax
   _REAL_ h, gox, goy, goz
   integer xm, ym, zm, xmymzm, level, nfocus
   integer narcdot, maxarc
   integer nbnd, nbndx, nbndy, nbndz
   integer outflag(natom)
   _REAL_ gcrd(3,natom), acrd(3,natom), radi(natom), radip3(natom)
   integer marc(*), m2narc(maxarc,*), fstarc(*), arcatm(2,*), dotarc(*)
   _REAL_ arccrd(3,*), savarc(3,*)
   integer atmsas((xm+2)*(ym+2)*(zm+2)), insas((xm+2)*(ym+2)*(zm+2))
   _REAL_ lvlset((xm+2)*(ym+2)*(zm+2)), mlvlset((xm+2)*(ym+2)*(zm+2)), zv((xm+2)*(ym+2)*(zm+2))
   _REAL_ epsx(xmymzm+ym*zm), epsy(xmymzm+xm*zm), epsz(xmymzm+xm*ym)
   integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer savbcopt(nfocus), saopt, sasopt

   ! Local variables

   integer ip, iatm, buf, nwarn, xmymzm_ext
   integer i, j, k
   integer newown
   integer ierr
   integer, allocatable :: atmx(:,:,:)
   integer, allocatable :: atmy(:,:,:)
   integer, allocatable :: atmz(:,:,:)
   integer, allocatable :: inmem(:,:,:)
   _REAL_ xi, yi, zi
   _REAL_ range0, range1, rh
   _REAL_, allocatable :: clsx(:,:,:)
   _REAL_, allocatable :: clsy(:,:,:)
   _REAL_, allocatable :: clsz(:,:,:)

   call pbtimer_start(PBTIME_PBEXMOL)

   ! initializations and local working arrays setup

   call pbtimer_start(PBTIME_PBEXMOL_SETUP)
    if ( sasopt == 0 ) then
      allocate(atmx(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      allocate(atmy(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      allocate(atmz(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      allocate(clsx(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      allocate(clsy(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      allocate(clsz(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      clsx = 1.d0; clsy = 1.d0; clsz = 1.d0 
      atmx = 0; atmy = 0; atmz = 0 
   end if

   if ( memopt > 0 ) then
      allocate(inmem(0:xm+1,0:ym+1,0:zm+1), stat = ierr )
      inmem = 0
   end if 

   xmymzm_ext = (xm+2)*(ym+2)*(zm+2)
   if (ipb /= 4 .and. ipb /= 5) then
      epsx(1:xmymzm+ym*zm) = epsout; epsy(1:xmymzm+xm*zm) = epsout; epsz(1:xmymzm+xm*ym) = epsout
   end if

   lvlset(1:xmymzm_ext) = 0.0d0; mlvlset(1:xmymzm_ext) = 0.0d0
   call pbtimer_stop(PBTIME_PBEXMOL_SETUP)

   ! in uniform dielectric systems, nothing to do here
   ! this is most likely for some reference state calculation or development work

   if ( epsin == epsout ) return

   ! in heterogeneous dielectrics systems ...

   ! labeling notation by insas():
   !
   ! -4: bulk solvent and salt
   ! -3: Stern layer
   ! -2: SAS layer + 4
   !
   ! if sasopt == 1, VDW/SAS surface
   ! +2: inside VDW/SAS
   !
   ! if sasopt == 0, SES surface
   ! -1: reentry region outside SES
   ! +1: reentry region inside SES
   ! +2: inside VDW

   newown = 1
   rh = 1.0d0/h

   ! CQ: turn on SAS surface when sasopt == 1
   if ( sasopt ==  1 ) then

      ! part a. reset atmsas and insas
      call pbtimer_start(PBTIME_PBEXMOL_PARTA)
      insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0
      call pbtimer_stop(PBTIME_PBEXMOL_PARTA)

      ! part b. mark grid points just a bit larger than SAS as -2

      ! WJ:
      ! 1. if use the level set-based SES, we first need to take care of the possibly very small
      ! solvent probe since we need level set function values on not just
      ! boundary grid points but its neighobrs as well. Adding 6 grids beyound
      ! dprob makes it safer for the later search of neighbors
      ! 2. this is also done for classical SES for comparison only ...

      call pbtimer_start(PBTIME_PBEXMOL_PARTB)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+6.0d0; xi = gcrd(1,iatm);yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTB)

      ! part c. mark grid points within SAS as 2, outside remains to be -2 from part b

      call pbtimer_start(PBTIME_PBEXMOL_PARTC)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( 2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTC)

   ! XP: turn on revised density-based surface when sasopt == 2
   ! WMBS: added membrane support

   else if ( sasopt == 2 ) then
      insas(1:xmymzm_ext) = -2
      call density_lvlset(memopt,outlvlset,outmlvlset,natom,xm,ym,zm,gox,goy,goz,h,dprob*rh,&
              cutfd,mzmin,mzmax,gcrd,radip3,poretype,poreradi,insas,inmem,lvlset,mlvlset)

   ! turn on SES for the finest grid when sasopt == 0
   ! Mengjuei: I've tested multiblock with ses enabled for the coarse grid. The
   ! accuracy was not improved by much. Not worth it.

   else if ( sasopt == 0 .and. level == nfocus ) then

      ! part a. reset atmsas and insas

      call pbtimer_start(PBTIME_PBEXMOL_PARTA)
      insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0
      call pbtimer_stop(PBTIME_PBEXMOL_PARTA)

      ! part b. mark grid points just a bit larger than SAS as -2

      ! WJ:
      ! 1. if use the level set-based SES, we first need to take care of the possibly very small
      ! solvent probe since we need level set function values on not just
      ! boundary grid points but its neighobrs as well. Adding 4 grids beyound
      ! dprob makes it safer for the later search of neighbors
      ! 2. this is also done for classical SES for comparison only ...

      call pbtimer_start(PBTIME_PBEXMOL_PARTB)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+4.0d0; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTB)

      ! part c. mark grid points within SAS as 1, outside remains to be -2 from part b

      call pbtimer_start(PBTIME_PBEXMOL_PARTC)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( 1, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTC)
    
      ! part d. mark grid points within VDW as 2, outside is 1 from part c
    
      call pbtimer_start(PBTIME_PBEXMOL_PARTD)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radi(iatm)
         if ( range1 == 0.0d0 ) cycle
         range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTD)
    
      ! part e. mark grid points in the contact region accessible to the solvent probe as -2
    
      call pbtimer_start(PBTIME_PBEXMOL_PARTE)
      call contact( insas, atmsas )
      call pbtimer_stop(PBTIME_PBEXMOL_PARTE)
    
      ! part f. mark grid points in the reentry region accessible to the solvent probe as -1
         
      call pbtimer_start(PBTIME_PBEXMOL_PARTF)
      buf = 2; range1 = dprob*rh + buf
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = narcdot, 1, -1
         xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
         call exresph( -1, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTF)

   ! when sasopt == 0 use modified VDW to approximate
   ! SES for the coarse grid
 
   else if ( sasopt == 0 .and. level /= nfocus ) then
       
      ! part a. reset atmsas and insas

      call pbtimer_start(PBTIME_PBEXMOL_PARTA)
      insas(1:xmymzm_ext) = -4; atmsas(1:xmymzm_ext) = 0
      call pbtimer_stop(PBTIME_PBEXMOL_PARTA)

      ! part b. mark grid points just a bit larger than VDW as -2

      call pbtimer_start(PBTIME_PBEXMOL_PARTB)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = radip3(iatm)
         if (range1 == 0.0d0 ) cycle
         range1 = (range1+dprob)*rh+4.0d0; xi = gcrd(1,iatm);yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exsasph( -2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTB)

      ! part c. mark volume within VDW as 2
      ! 1. since outside is -2, there are only contact fractional edges
      ! 2. we are using the modified VDW radii that have been
      ! augmented by radinc in sa_driver()

      call pbtimer_start(PBTIME_PBEXMOL_PARTD)
      zv(1:xmymzm_ext) = 9999.0d0
      do iatm = 1, natom
         range1 = radip3(iatm)
         if ( range1 == 0.0d0 ) cycle
         if (ifcap == 5 .and. outflag(iatm) == 1) cycle
         range1 = range1*rh; xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
         call exvwsph( 2, insas, atmsas, zv )
      end do
      call pbtimer_stop(PBTIME_PBEXMOL_PARTD)
 
   end if

   ! Here are the extras not needed by every option:

   ! part a.
   ! For classical PB surface, ipb == 1 or 2, the membrane is modeled as a slab
   ! for sasopt = 0 and 1, i.e. SES and SAS.

   if ( memopt > 0 .and. sasopt /= 2 .and. (ipb == 1 .or. ipb == 2) ) then
      ! step 1. Mark solvent grid points as membrane if they are within the
      ! membrane zone, i.e. with mzmin < z < mzmax
      call exvwslab( 1,insas,inmem )

      ! step 2. Find the pore that goes through the membrane by default
      ! unless this is not turned on.
      if ( poretype /= 0 ) call findpore( insas,inmem )
   end if

   ! part b. 
   ! identify irregular grid points. Also identify atom/probe owner if a
   ! geometry based method (sasopt = 0/1) is used.
   !
   ! The eps boundary grid location routine only works for globular proteins.
   ! The physics is unclear on how to do it for a slab in a periodic box.
   !
   ! save boundary edges for eneopt=2 energy and forces
   ! CQ: always call epsbnd because bcopt == 6 requires this 
   ! WJ: epsbnd has been revised. should be free of warning in normal calls

   call pbtimer_start(PBTIME_PBEPSBND)
   if ( memopt == 0 ) then
      call epsbnd( atmsas, insas )
   end if

   ! part c.
   ! WJ: if requested, set up the level set function as signed distance to SES
   ! RL: sasopt == 1 is now implemented for both SAS and VDW (i.e. SAS/dprob=0)

   if ( ( ipb == 2 .or. ipb == 4 .or. ipb == 5 ) .and. sasopt /= 2 ) then
      lvlset(1:xmymzm_ext) = -dble(sign(9999, insas(1:xmymzm_ext)))
      if      ( sasopt == 0 ) then
         if ( level == nfocus ) call assignlvlset( natom,radi(1:natom),atmsas,insas,lvlset )
         if ( level /= nfocus ) call assignlvlset( natom,radip3(1:natom),atmsas,insas,lvlset )
      else if ( sasopt == 1 ) then
         call assignlvlset( natom,radi(1:natom)+dprob,atmsas,insas,lvlset )
      else
         write(6,'(a)') 'PB Bomb in pb_exmol: Incompatible option for lvlset assignment'
         call mexit(6,1)
      end if
   end if

   if ( (ipb == 2 .or. ipb == 4 .or. ipb == 5) .and. memopt == 0 .and. level == nfocus ) then 
      call irreg_init(xm,ym,zm,nbnd,h,gox,goy,goz,insas,iepsav,lvlset)
      call irreg_grd(xm,ym,zm,nbnd,h,iepsav,lvlset)
      if ( sasopt == 2 ) &
      call bndatm(nbnd,natom,xm,ym,zm,h,dprob*rh,cutfd,gcrd,radip3,iepsav,lvlset)
   end if
   call pbtimer_stop(PBTIME_PBEPSBND)

   ! part d.
   ! set up epsx, epsy, and epsz maps. boundary edges for eneopt=1 and
   ! frcopt=1&3 are saved inside. this is not necessary for the newer solvers
   ! that enforce interface conditions
   ! RL: both atomic and membrane dielectric grid edges are handled here

   call pbtimer_start(PBTIME_PBEPSMAP)
   if ( ipb == 1 .or. ipb == 2 ) then
      call epsmap( ipb,memopt,insas,atmsas,inmem,lvlset,epsx,epsy,epsz )
   end if
   call pbtimer_stop(PBTIME_PBEPSMAP)

   ! part e.
   ! CQ: calculate SAS or SES surface area based on the finite-difference data
   ! structure if requested

   call pbtimer_start(PBTIME_PBCALSA)
   if ( level == nfocus .and. sasopt <= 1 .and. (.not. ligand .or. .not. multiblock)  ) then
      if ( abs(saopt) == 1 ) &
         call calc_sa1(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                 iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt)
      if ( abs(saopt) == 2 ) &
         call calc_sa2(acrd,xm,ym,zm,xmymzm,nbnd,iepsav,iepsavx,iepsavy,&
                 iepsavz,gox,goy,goz,h,sasopt,smoothopt, &
                 epsin,epsout,insas,epsx,epsy,epsz)
   end if
   call pbtimer_stop(PBTIME_PBCALSA)

   ! part f.
   ! cleaning up ...

   call pbtimer_start(PBTIME_PBEXMOL_SETUP)
   if ( sasopt == 0 ) then
      deallocate(atmx, stat = ierr )
      deallocate(atmy, stat = ierr )
      deallocate(atmz, stat = ierr )
      deallocate(clsx, stat = ierr )
      deallocate(clsy, stat = ierr )
      deallocate(clsz, stat = ierr )
   end if
   if ( memopt > 0 ) then
      deallocate(inmem, stat = ierr )
   end if 
   call pbtimer_stop(PBTIME_PBEXMOL_SETUP)

   call pbtimer_stop(PBTIME_PBEXMOL)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic sasurf
subroutine exsasph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom  (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
 
   ! Passed variables

   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d, d2, r

   r = rh*radi(iatm)

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2; d = sqrt(d2)
               if ( d > r ) then
                  d = d - r
                  if ( insph(i,j,k) == dielsph ) then
                     if ( d < dst(i,j,k) ) then
                        inatm(i,j,k) = iatm; dst(i,j,k) = d
                     end if
                     cycle
                  end if
                  insph(i,j,k) = dielsph;
                  inatm(i,j,k) = iatm; dst(i,j,k) = d
               else
                  if ( insph(i,j,k) == dielsph ) cycle
                  insph(i,j,k) = dielsph; inatm(i,j,k) = iatm; dst(i,j,k) = 0.0d0
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk

end subroutine exsasph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none
 
   ! Passed variables
    
   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3
   _REAL_ projection(1:3),tmp(1:3),xl(1:3),dist,d2
   integer iatml

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               if ( insph(i,j,k) == dielsph ) then
                  if ( level /= nfocus ) cycle
                  iatml = inatm(i,j,k)
                  xl = gcrd(1:3,iatml)
                  tmp(1) = i - xl(1)
                  tmp(2) = j - xl(2)
                  tmp(3) = k - xl(3)
                  dist = sqrt(sum(tmp*tmp))
                  if ( dist == 0.d0 ) then
                     inatm(i,j,k) = iatm
                     cycle
                  end if
                  tmp = tmp / dist
                  projection = xl+tmp*radi(iatml)*rh
                  tmp(1) = projection(1)-xi
                  tmp(2) = projection(2)-yi 
                  tmp(3) = projection(3)-zi
                  d2 = sum(tmp*tmp)
                  if ( d2 < range1*range1 ) then
                     inatm(i,j,k) = iatm
                  end if
               else
                  insph(i,j,k) = dielsph
                  inatm(i,j,k) = iatm
               end if
            end do  ! i = lowi, highi
             
         end if  ! ( range3 >= 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
   
end subroutine exvwsph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark contact grid points between vdw and sas surfaces
subroutine contact( insas,atmsas )
    
   implicit none
 
   integer insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   integer i, j, k, iatm, ii, jj, iarc, inside
   _REAL_ xg(3), xi(3), xj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji

   _REAL_, parameter :: small = 0.01d0
    
   do k = 0, zm+1; do j = 0, ym+1; do i = 0, xm+1
       
      if ( insas(i,j,k) < 0 ) cycle
       
      xg(1) = gox + i*h; xg(2) = goy + j*h; xg(3) = goz + k*h
       
      ! this is the atom that marked this grid within sasrf, so it will be the
      ! grid's contact atom if it is marked so.
       
      iatm = atmsas(i,j,k)
      xi(1) = acrd(1,iatm); xi(2) = acrd(2,iatm); xi(3) = acrd(3,iatm)
       
      ! go through all arcs that this atom generates in circle()

      if ( fstarc(iatm) == 0 ) then
         write(6,'(a)') 'PB Bomb in contact(): fstarc=0'
         write(6,'(a,i5,3f10.3)') 'atom info:', iatm, xi(1:3), radi(iatm), radi(iatm) + dprob
         write(6,'(a,f10.3,4i5)') 'grid info:', xg(1:3), i, j, k, insas(i,j,k) 
         call mexit(6,1)
      end if
      inside = -2
      do ip = 1, marc(iatm)
         iarc = m2narc(ip,iatm)

         ! generated by outer loop, i.e. the atom is iatm in circle()
          
         if ( iarc >= fstarc(iatm) ) then
            jj = arcatm(1,iarc)
            xj(1) = acrd(1,jj)
            xj(2) = acrd(2,jj)
            xj(3) = acrd(3,jj)
            cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
          
         ! generated by inner loop, i.e. the atom is jatm in circle()
          
         else
            xj = xi
            ii = arcatm(2,iarc)
            xj(1) = acrd(1,ii)
            xj(2) = acrd(2,ii)
            xj(3) = acrd(3,ii)
            cosaji = savarc(1,iarc); cosaij = savarc(2,iarc)
         end if
         rxij = savarc(3,iarc)
         xij = rxij*(xj - xi)
          
         dx = xg - xi; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         dx = xg - xj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))
          
         ! if gij < aij .and. gji < aji, this is a reentry grid
          
         if ( cosgij <= cosaij .or. cosgji <= cosaji ) cycle
         if ( insas(i,j,k) == 1 ) then
            inside = 1
         else if ( insas(i,j,k) == 2 ) then
            inside = 1
         end if
         exit
      end do
       
      if ( inside == -2 .and. insas(i,j,k) == 2 ) inside = 2
      insas(i,j,k) = inside
       
   end do; end do; end do
   
end subroutine contact
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within reentry surf
subroutine exresph( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a renentry sphere (dielectric constant dielsph)
   ! of index iatm as dielsph. Modified from UHBD (Comp. Phys. Comm. 91:57-95,
   ! 1995) routines excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)

   ! Local variables

   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2, front, aa

   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk
       
      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj
          
         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then
             
            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi

               if ( insph(i,j,k) == 2 ) cycle
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == -2 ) then
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) == 1 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) == 1 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) == 1 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
               else if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) - 1.d-9 ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
               else if ( insph(i,j,k) == 1 ) then
                  if ( d2 < (range1 - buf)**2 ) then
                     insph(i,j,k) = dielsph;
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  if ( newown == 1 ) then
                  if ( i > 0 ) then
                     if ( insph(i-1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi-sqrt(front)-REAL(i-1)
                        if ( aa >= 0.d0 .and. aa < clsx(i-1,j,k) ) then
                           clsx(i-1,j,k) = aa
                           atmx(i-1,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( i < xm+1 ) then
                     if ( insph(i+1,j,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = xi+sqrt(front)-REAL(i)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsx(i,j,k) ) then
                           clsx(i,j,k) = 1.d0-aa
                           atmx(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j > 0 ) then
                     if ( insph(i,j-1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi-sqrt(front)-REAL(j-1)
                        if ( aa >= 0.d0 .and. aa < clsy(i,j-1,k) ) then
                           clsy(i,j-1,k) = aa
                           atmy(i,j-1,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( j < ym+1 ) then
                     if ( insph(i,j+1,k) > 0 ) then
                     front = (range1-buf)**2-(zi-k)**2-(xi-i)**2 
                     if ( front >= 0.d0 ) then
                        aa = yi+sqrt(front)-REAL(j)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsy(i,j,k) ) then
                           clsy(i,j,k) = 1.d0-aa
                           atmy(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k > 0 ) then
                     if ( insph(i,j,k-1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi-sqrt(front)-REAL(k-1)
                        if ( aa >= 0.d0 .and. aa < clsz(i,j,k-1) ) then
                           clsz(i,j,k-1) = aa
                           atmz(i,j,k-1) = iatm
                        end if
                     end if
                     end if
                  end if
                  if ( k < zm+1 ) then
                     if ( insph(i,j,k+1) > 0 ) then
                     front = (range1-buf)**2-(xi-i)**2-(yi-j)**2 
                     if ( front >= 0.d0 ) then
                        aa = zi+sqrt(front)-REAL(k)
                        if ( aa <= 1.d0 .and. 1.d0-aa < clsz(i,j,k) ) then
                           clsz(i,j,k) = 1.d0-aa
                           atmz(i,j,k) = iatm
                        end if
                     end if
                     end if
                  end if
                  end if 
                  else if ( d2 < dst(i,j,k) - 1.d-9 ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
               end if
!              insph(i,j,k) = dielsph;
!              inatm(i,j,k) = iatm; dst(i,j,k) = d2

            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
   
end subroutine exresph
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within an rectangular slab in parallel to the xy plane
subroutine exvwslab( dielslab,insas,inslab )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within a slab (dielectric constant dielsph) between
   ! zmin=range0 and zmax=range1 (in the FD grid unit and frame) as dielpsh if
   ! they are outside the solute molecules.
   !
   ! Note the difference with the routine for atoms, which also saves atom ID for 
   ! downstream force calculations.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   implicit none
 
   ! Passed variables

   integer dielslab
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   integer inslab(0:xm+1,0:ym+1,0:zm+1)
 
   ! Local variables
    
   integer lowk, highk

   lowk = max(0,ceiling(mzmin)); highk = min(zm+1,floor(mzmax))
   do k = lowk, highk
      do j = 0, ym+1; do i = 0, xm+1
         if ( insas(i,j,k) < 0 ) inslab(i,j,k) = 1
         !if ( insas(i,j,k) < 0 ) insas(i,j,k) = 0 ! this is for exporting to the dx file
      end do; end do
   end do
   !call gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,insas(1:xm,1:ym,1:zm),'insas.dx',100,'insas')

end subroutine exvwslab             
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the pore region that is accessible to water
subroutine findpore( insas,inslab )
   !
   ! 
   !

   implicit none
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   integer inslab(0:xm+1,0:ym+1,0:zm+1)

   ! local varaibles

   integer nzone, ier
   integer i,j,k,i0,j0,k0
   integer lowk, highk

   integer, allocatable :: kzone(:,:,:)
  
   if ( allocated(kzone) ) then
      deallocate(kzone, stat = ier); REQUIRE(ier==0)
   end if
   allocate(kzone(0:xm+1,0:ym+1,0:zm+1), stat = ier); REQUIRE(ier==0)

   kzone = -1 ! flags showing it's not visited, otherwise, it's the zone label.

   ! solute interior is zone # 0 whether it's connected or not
   ! find solvent accessible membrane zone in the system, which is zone # 1
   ! assuming we start from the grid edge, which belongs to zone # 1
    
   nzone = 0 ! no. of disconnected zones
   lowk = max(0,ceiling(mzmin)); highk = min(zm+1,floor(mzmax))
   do k = lowk, highk
      do j = 0, ym+1
      do i = 0, xm+1
         if ( kzone(i,j,k) /= -1 ) cycle
         if ( insas(i,j,k) < 0 ) then
            nzone = nzone + 1
            kzone(i,j,k) = nzone
            i0 = i; j0 = j; k0 = k 
            call walk(0,i0,j0,k0,xm,ym,zm,insas,kzone,nzone,lowk,highk)
         else
            kzone(i,j,k) = 0
         end if
      end do
      end do
   end do

   ! now revise the membrane lables so that the grid points that are neither
   ! zone 0 (solute) or zone 1 (membrane) are reset back to solvent.
    
   do k = lowk, highk
      do j = 0, ym+1; do i = 0, xm+1
         if ( kzone(i,j,k) > 1 ) inslab(i,j,k) = 0
         !if ( insas(i,j,k) < 0 ) insas(i,j,k) = 0 ! this is for exporting to the dx file
      end do; end do
   end do
   !call gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,inslab(1:xm,1:ym,1:zm),'inmem.dx',101,'inmem')
   !call gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,insas(1:xm,1:ym,1:zm),'insas.dx',102,'insas')

end subroutine findpore
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ the depth first search algorithm
recursive subroutine walk(al,i,j,k,xm,ym,zm,phi,kzone,nzone,lowk,highk)

   implicit none
   integer al   
   integer i,j,k
   integer xm,ym,zm
   integer nzone
   integer lowk,highk
   integer phi(0:xm+1,0:ym+1,0:zm+1)
   integer kzone(0:xm+1,0:ym+1,0:zm+1)

   kzone(i,j,k) = nzone

   if ( i<xm+1  .and. kzone(i+1,j,k) == -1 .and. phi(i+1,j,k) < al ) &
      call walk(al,i+1,j,k,xm,ym,zm,phi,kzone,nzone,lowk,highk)
   if ( i>0     .and. kzone(i-1,j,k) == -1 .and. phi(i-1,j,k) < al ) &
      call walk(al,i-1,j,k,xm,ym,zm,phi,kzone,nzone,lowk,highk)
   if ( j<ym+1  .and. kzone(i,j+1,k) == -1 .and. phi(i,j+1,k) < al ) &
      call walk(al,i,j+1,k,xm,ym,zm,phi,kzone,nzone,lowk,highk)
   if ( j>0     .and. kzone(i,j-1,k) == -1 .and. phi(i,j-1,k) < al ) &
      call walk(al,i,j-1,k,xm,ym,zm,phi,kzone,nzone,lowk,highk)
   if ( k<highk .and. kzone(i,j,k+1) == -1 .and. phi(i,j,k+1) < al ) &
      call walk(al,i,j,k+1,xm,ym,zm,phi,kzone,nzone,lowk,highk)
   if ( k>lowk  .and. kzone(i,j,k-1) == -1 .and. phi(i,j,k-1) < al ) &
      call walk(al,i,j,k-1,xm,ym,zm,phi,kzone,nzone,lowk,highk)

end subroutine walk
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd ( atmsas,insas )
    
   use poisson_boltzmann, only: savgox, savgoy, savgoz
   implicit none
    
   integer atmsas(0:xm+1,0:ym+1,0:zm+1), insas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   logical boundary
   integer i, j, k, clstmp
   integer iarc

   nwarn = 0
   nbnd = 0
   do k = 1, zm; do j = 1, ym; do i = 1, xm

      ! set up condition for a boundary grid point

      boundary = .false.
      if ( (insas(i,j,k) ==  1 .or. insas(i,j,k) ==  2) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) == -1 .or. insas(i,j,k) == -2) .and.&
           (insas(i-1,j,k) ==  1 .or. insas(i-1,j,k) ==  2 .or. insas(i+1,j,k) ==  1 .or.&
            insas(i+1,j,k) ==  2 .or. insas(i,j-1,k) ==  1 .or. insas(i,j-1,k) ==  2 .or.&
            insas(i,j+1,k) ==  1 .or. insas(i,j+1,k) ==  2 .or. insas(i,j,k-1) ==  1 .or.&
            insas(i,j,k-1) ==  2 .or. insas(i,j,k+1) ==  1 .or. insas(i,j,k+1) ==  2 .or.&
            insas(i-1,j,k) ==  0 .or. insas(i+1,j,k) ==  0 .or. insas(i,j-1,k) ==  0 .or.&
            insas(i,j+1,k) ==  0 .or. insas(i,j,k-1) ==  0 .or. insas(i,j,k+1) ==  0) ) then
            boundary = .true.
      else if ( (insas(i,j,k) ==  0) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2 &
           ) ) then 
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k

      if ( ifcap /= 0 .and. ifcap /= 5 ) then
         iepsav(4,nbnd) = 0
         cycle
      end if

      ! if we are using a geometry based surface, i.e. SES or SAS/VDW.
 
      if ( sasopt < 2 ) then  

         ! for a grid point in contact region +/- 2 or in a solvent probe, simply use the atom/probe that
         ! marks it

         clstmp = 0
         if ( abs(insas(i,j,k)) == 2 ) then
            clstmp = atmsas(i,j,k)
            if ( clstmp == 0 ) then
               if ( pbverbose .and. level == nfocus ) then
                  write(6, '(a,4i4)') &
                  'PB Warning in epsbnd(): No neighbor found for exposed boundary grid',&
                  i, j, k, insas(i,j,k)
               end if
               nwarn = nwarn + 1
            end if
         else if ( insas(i,j,k) == -1 ) then
            clstmp = -atmsas(i,j,k)
            if ( clstmp == 0 ) then
               if ( pbverbose .and. level == nfocus ) then 
                  write(6, '(a,4i4)') &
                  'PB Warning in epsbnd(): No neighbor found for exposed boundary grid',&
                   i, j, k, insas(i,j,k)
               end if
               nwarn = nwarn + 1
            end if
 
         ! for a buried reentry grid point, find the atom that marked its neighoring exposed reentry
         ! grid points. Note that this may not be possible when grid spacing is large
 
         else if ( insas(i,j,k) == 1 ) then
            clstmp = -atmsas(i,j,k)
            if ( clstmp == 0 ) then
               write(6,'(a)') 'PB Info: Close contact cannot be found. fndcls() is called'
               clstmp = fndcls( i, j, k, insas, atmsas )
               nwarn = nwarn + 1
            end if

         end if
         iepsav(4,nbnd) = clstmp

      end if

   end do; end do; end do
   if ( nwarn > 0 ) then
      if ( pbverbose .and. level == nfocus ) write(6, '(a,i4)') &
      'PB Warning in epsbnd(): No neighbor found for boundary grids total:', nwarn
   end if
 
end subroutine epsbnd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls( i,j,k,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls, i, j, k
   integer  insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - dprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls = clsatm(ip)
   end do

end function fndcls
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign level set function for grid points nearby the interface
subroutine assignlvlset( natom,radi,atmsas,insas,u )

   implicit none

   ! passed variables

   integer natom
   _REAL_ radi(natom)
   integer atmsas(0:xm+1,0:ym+1,0:zm+1), insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer i,j,k,l,buffer
   _REAL_ dist,d2,d

   ! the basic idea is to handle grid points nearby the surface. Apparently the
   ! way it is set up repeats assignment for most visited grid points, so an if
   ! is inserted to skip those done grid points.
 
   rh = 1.0d0 / h
   dist = dprob*rh
   buffer = 2 
   do l = 1, nbnd
      do k = iepsav(3,l) - buffer, iepsav(3,l) + buffer
         do j = iepsav(2,l) - buffer, iepsav(2,l) + buffer
            do i = iepsav(1,l) - buffer, iepsav(1,l) + buffer
               if ( abs(u(i,j,k)) /= 9999.0d0 ) cycle ! this is done already
               if ( atmsas(i,j,k) == 0 ) then
                  write(6,'(a,4i5)') 'PB Bomb in assignlvlset(): no atmsas', i,j,k, insas(i,j,k)
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l))
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l)+1,iepsav(2,l),iepsav(3,l))
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l)-1,iepsav(2,l),iepsav(3,l))
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)+1,iepsav(3,l))
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l)-1,iepsav(3,l))
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)+1)
                  write(6,'(4i6)') iepsav(1:3,l),insas(iepsav(1,l),iepsav(2,l),iepsav(3,l)-1)
                  call mexit(6,1)
               end if
               if ( abs(insas(i,j,k)) == 2 ) then 
                  iatm = atmsas(i,j,k)
                  xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d = sqrt(d2)
                  u(i,j,k) = - radi(iatm)*rh + d 
               else if ( abs(insas(i,j,k)) == 1 ) then
                  iatm = atmsas(i,j,k)
                  xi = (arccrd(1,iatm) - gox)*rh; yi = (arccrd(2,iatm) - goy)*rh; zi = (arccrd(3,iatm) - goz)*rh
                  d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2 ; d =sqrt(d2)
                  u(i,j,k) = dist - d
               else
                  write(6,'(a,4i5)') 'PB Bomb in assignlvlset(): illegal insas flag', i,j,k, insas(i,j,k)
                  call mexit(6,1)
               end if
            end do
         end do
      end do
   end do

   ! dx based visualization of the surface
   !call gen_dx_file(xm,ym,zm,h,gox,goy,goz,u(1:xm,1:ym,1:zm),'lvlset.dx',100,'lvlset') 

end subroutine assignlvlset
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap( ipb,memopt,insas,atmsas,inmem,u,epsx,epsy,epsz )

   implicit none

   ! passed variables

   integer ipb, memopt
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   integer atmsas(0:xm+1,0:ym+1,0:zm+1)
   integer inmem(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsx(0:xm,1:ym,1:zm)
   _REAL_ epsy(1:xm,0:ym,1:zm)
   _REAL_ epsz(1:xm,1:ym,0:zm)

   ! local variables

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   integer aa, bb, cc, dd
   integer x_flag, y_flag, z_flag
   _REAL_ epsint, epsint0
   _REAL_ frac
   _REAL_ epslow, epshigh

   if ( sasopt == 2 .and. memopt > 0 ) then
      write(6,'(a)') ' PB Bomb in pb_epsmap(): The density function option for'
      write(6,'(a)') ' the membrane is not supported in this release.'
      write(6,'(a)') ' Please use the SES option instead, sasopt=0.'
      call mexit(6,1)
   end if

   ! set default value for simple harmonic average and the coarse grid level

   epsint0 = 2.0d0*epsin*epsout/(epsin+epsout)
   epsint = epsint0

   ! initialize fraction edge counters

   nbndx = 0; nbndy = 0; nbndz = 0
   x_flag = 0; y_flag = 0; z_flag = 0 ! local copies for the level set version

   do k = 0, zm; do j = 0, ym; do i = 0, xm

      ! x-edges:

      a  = insas(i,j,k);  b  = insas(i+1,j,k)
      a1 = atmsas(i,j,k); b1 = atmsas(i+1,j,k)
      if ( a == 0 .or. b == 0 ) go to 10 ! temperary skip for sasopt=2

      aa = 0; bb = 0
      if ( memopt > 0 .and. sasopt == 0 ) then
         aa = inmem(i,j,k); bb = inmem(i+1,j,k)
      end if

      ! I. for grid edges acrossing the membrane boundaries
      !
      ! a. edges flanked by two "1" points are epsmem edges
      ! b. edges flanked by one "1" point and one positive point are fractional solute edges 
      ! c. edges flanked by one "1" point and one negative point are fractional solvent edges 

      if ( aa == 1 .or. bb == 1 ) then

         if ( j == 0 .or. k == 0 ) then ! note this edge is defined
            continue
         else if ( aa == 1 .and. bb == 1 ) then ! case a
            ! case a
            ! this edge has nothing to do with either solute or solvent
            epsx(i,j,k) = epsmem
         else if ( a > 0 .or. b > 0 ) then ! case b
            ! case b
            ! the membrane is treated as a solvent phase by setting epsout = epsmem.
            call epsfracx(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsmem)
            epsx(i,j,k) = epsint
         else if ( a < 0 .or. b < 0 ) then ! case c
            ! case c
            ! since the membrane is a slab parallel to the xy plane, there can't be
            ! such an fractional edge.
            continue
         else
            ! faulty combo
            write(6,'(a)') 'PB Bomb in pb_epsmap(): inconsistent grid labels for x membrane edges'
            call mexit(6,1)
         end if

      ! II. for rest of the grid eges, the single-solvent solvated situation
      ! applies.
      !
      ! a. edges flanked by two positive points are epsin edges
      ! b. edges flanked by two negative points are epsout edges and
      ! don't have to be reset
      ! c. edges flanked by two points of different signs are fraction edges

      else

         if ( j == 0 .or. k == 0 ) then ! note this edge is not declared
            continue
         else if ( sign(a,b) == a ) then
            if ( a > 0 ) then
               epsx(i,j,k) = epsin
            end if
         else
            if ( level == nfocus ) then
               ! for all smoothopt's do the edge fraction calculation
               if ( (ipb == 1 .or. ipb == 2 .or.&
                  (ipb == 3 .and. (a == 2 .or. b == 2)) ) .and. sasopt /= 2 ) & 
                  call epsfracx  (i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout  )
               if ( ((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. b == 1)) ) .and. &
                  nbndx > x_flag) .or. (ipb == 2  .and. sasopt == 2) ) then
                  call epsfracx_r(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout,u)
                  x_flag = x_flag + 1
               end if
               ! this option is meant to reproduce the delphi behavior
               if ( smoothopt == 2 ) then
                  if ( epsint > epsint0 ) then
                     epsint = epsout
                  else
                     epsint = epsin
                  end if
               end if
               ! this dump option is for testing only
               if ( smoothopt == 0 ) then 
                  epsint = epsint0
               end if
            end if
            epsx(i,j,k) = epsint
         end if

      end if

      ! y-edges:

10    continue
      c  = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( a == 0 .or. c == 0 ) go to 20 ! temperary skip for sasopt=2

      cc = 0
      if ( memopt > 0 .and. sasopt == 0 ) cc = inmem(i,j+1,k)

      ! I. for grid edges acrossing the membrane boundaries
      !
      ! a. edges flanked by two "1" points are epsmem edges
      ! b. edges flanked by one "1" point and one positive point are fractional solute edges 
      ! c. edges flanked by one "1" point and one negative point are fractional solvent edges 

      if ( aa == 1 .or. cc == 1 ) then

         if ( i == 0 .or. k == 0 ) then ! note this edge is not declared
            continue
         else if ( aa == 1 .and. cc == 1 ) then
            ! case a
            ! this edge has nothing to do with either solute or solvent
            epsy(i,j,k) = epsmem
         else if ( a > 0 .or. c > 0 ) then
            ! case b
            ! the membrane is treated as a solvent phase by setting epsout = epsmem.
            call epsfracy(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsmem)
            epsy(i,j,k) = epsint
         else if ( a < 0 .or. c < 0 ) then
            ! case c
            ! since the membrane is a slab parallel to the xy plane, there can't be
            ! such an fractional edge.
            continue
         else
            ! faulty combo
            write(6,'(a)') 'PB Bomb in pb_epsmap(): inconsistent grid labels for y membrane edges'
            call mexit(6,1)
         end if

      ! II. for rest of the grid eges, the single solvent solvated situation
      ! applies.
      !
      ! a. edges flanked by two positive points are epsin edges
      ! b. edges flanked by two negative points are epsout edges and
      ! don't have to be reset
      ! c. edges flanked by two points of different signs are fraction edges

      else

         if ( i == 0 .or. k == 0 ) then
            continue
         else if ( sign(a,c) == a ) then
            if ( a > 0 ) then
               epsy(i,j,k) = epsin
            end if
         else
            if ( level == nfocus ) then
               ! for all smoothopt's do the edge fraction calculation
               if ( (ipb == 1 .or. ipb == 2 .or.&
                    (ipb == 3 .and. (a == 2 .or. c == 2)) ) .and. sasopt /= 2 ) &
                  call epsfracy  (i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout  )
               if ( ((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. c == 1)) ) .and.&
                    nbndy > y_flag) .or. (ipb == 2 .and. sasopt == 2) ) then
                  call epsfracy_r(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout,u)
                  y_flag = y_flag + 1
               end if
               ! this option is meant to reproduce the delphi behavior
               if ( smoothopt == 2 ) then
                  if ( epsint > epsint0 ) then
                     epsint = epsout
                  else
                     epsint = epsin
                  end if
               end if
               ! this dump option is for testing only
               if ( smoothopt == 0 ) then 
                  epsint = epsint0
               end if
            end if
            epsy(i,j,k) = epsint
         end if

      end if

      ! z-edges:

20    continue
      d  = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( a == 0 .or. d == 0 ) go to 30 ! temperary skip for sasopt=2 

      dd = 0
      if ( memopt > 0 .and. sasopt == 0 ) dd = inmem(i,j,k+1)

      ! I. for grid edges acrossing the membrane boundaries
      !
      ! a. edges flanked by two "1" points are epsmem edges
      ! b. edges flanked by one "1" point and one positive point are fractional solute edges 
      ! c. edges flanked by one "1" point and one negative point are fractional solvent edges 

      if ( aa == 1 .or. dd == 1 ) then

         if ( i == 0 .or. j == 0 ) then ! note this edge is not declared
            continue
         else if ( aa == 1 .and. dd == 1 ) then
            ! case a
            ! this edge has nothing to do with either solute or solvent
            epsz(i,j,k) = epsmem
         else if ( a > 0 .or. d > 0 ) then
            ! case b
            ! the membrane is treated as a solvent phase by setting epsout = epsmem.
            call epsfracz(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsmem)
            epsz(i,j,k) = epsint
         else
            ! case c
            ! the membrane is treated as a solute phase by setting epsin = epsmem.
            ! note too this is a brand-new situation that hasn't been taken care of
            ! by exising methods developed for single-solvent solution systems.
            ! However this is much simpler than packed spheres in molecules.
            if ( aa == 1 .and. d < 0 ) then
               frac = REAL(k+1) - mzmax
            else if ( dd == 1 .and. a < 0 ) then
               frac = mzmin - REAL(k)
            else
               ! faulty combo
               write(6,'(a)') 'PB Bomb in pb_epsmap(): inconsistent grid labels for z membrane edges'
               call mexit(6,1)
            end if
            epsint = (epsout*epsmem)/(epsmem*(1.0d0-frac) + epsout*frac)
            epsz(i,j,k) = epsint
         end if

      ! II. for rest of the grid eges, the single solvent solvated situation
      ! applies.
      !
      ! a. edges flanked by two positive points are epsin edges
      ! b. edges flanked by two negative points are epsout edges and
      ! don't have to be reset
      ! c. edges flanked by two points of different signs are fraction edges

      else

         if ( i == 0 .or. j == 0 ) then
            continue
         else if ( sign(a,d) == a ) then
            if ( a > 0 ) then
               epsz(i,j,k) = epsin
            end if
         else
            if ( level == nfocus ) then
               ! for all smoothopt's do the edge fraction calculation
               if ( (ipb == 1 .or. ipb == 2 .or.&
                    (ipb == 3 .and. (a == 2 .or. d == 2)) ) .and. sasopt /= 2 ) &
                  call epsfracz  (i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout  )
               if ( ((ipb == 2 .or. (ipb == 3 .and. ( a == 1 .or. d == 1)) ) .and. nbndz > z_flag) .or.&
                    (ipb == 2 .and. sasopt == 2) ) then
                  call epsfracz_r(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout,u)
                  z_flag = z_flag + 1
               end if
               ! this option is meant to reproduce the delphi behavior
               if ( smoothopt == 2 ) then
                  if ( epsint > epsint0 ) then
                     epsint = epsout
                  else
                     epsint = epsin
                  end if
               end if
               ! this dump option is for testing only
               if ( smoothopt == 0 ) then 
                  epsint = epsint0
               end if
            end if
            epsz(i,j,k) = epsint
         end if

      end if

30    continue
   end do; end do; end do

   ! checking the sanity of epsx, epsy, and epsz
   ! apparently min(epsin,epsmem,epsout) <= epsxyz <= max(epsin,epsmem,epsout)

   epslow = min(epsin,epsmem,epsout)*0.999d0
   epshigh = max(epsin,epsmem,epsout)*1.001d0
   do k = 0, zm; do j = 0, ym; do i = 0, xm

      if ( j == 0 .or. k == 0 ) then
         continue
      else if ( epsx(i,j,k) < epslow .or. epsx(i,j,k) > epshigh ) then
         write(6,'(a,3i10)') 'PB Bomb in epsmap(): epsx out of range', i,j,k
         write(6,'(a,3f12.4)') 'epsx: ',epsx(i,j,k)/min(epsin,epsmem,epsout)
         call mexit(6,1)
      end if
      if ( i == 0 .or. k == 0 ) then
         continue
      else if ( epsy(i,j,k) < epslow .or. epsy(i,j,k) > epshigh ) then
         write(6,'(a,3i10)') 'PB Bomb in epsmap(): epsy out of range', i,j,k
         write(6,'(a,3f12.4)') 'epsy:',epsin,epsout,epsy(i,j,k)/min(epsin,epsmem,epsout)
         call mexit(6,1)
      end if
      if ( i == 0 .or. j == 0 ) then
         continue
      else if ( epsz(i,j,k) < epslow .or. epsz(i,j,k) > epshigh ) then
         write(6,'(a,3i10)') 'PB Bomb in epsmap(): epsz out of range', i,j,k
         write(6,'(a,3f12.4)') 'epsz:',epsz(i,j,k)/min(epsin,epsmem,epsout)
         call mexit(6,1)
      end if
   end do; end do; end do

   ! here is the new dx debugging outputs
!  epsx = epsx/min(epsin,epsmem,epsout)
!  epsy = epsy/min(epsin,epsmem,epsout)
!  epsz = epsz/min(epsin,epsmem,epsout)
!  call gen_dx_file(xm,ym,zm,h,gox,goy,goz,epsx(1:xm,1:ym,1:zm),'epsx.dx',101,'epsx')
!  call gen_dx_file(xm,ym,zm,h,gox,goy,goz,epsy(1:xm,1:ym,1:zm),'epsy.dx',102,'epsy')
!  call gen_dx_file(xm,ym,zm,h,gox,goy,goz,epsz(1:xm,1:ym,1:zm),'epsz.dx',103,'epsz')
!  stop
end subroutine epsmap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = range3 + xi - REAL(i)
         fedgex(nbndx) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndx = nbndx + 1
         aa = range3 - xi + REAL(i+1)
         fedgex(nbndx) = 1.0d0 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + xi - REAL(i)
         xg(1) = gox + h*i + h*aa
         xg(2) = goy + h*j
         xg(3) = goz + h*k
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndx = nbndx + 1
            fedgex(nbndx) =  aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
            iepsavx(4,nbndx) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         nbndx = nbndx + 1
         aa = xi - range3 - REAL(i)
         fedgex(nbndx) = aa
         if ( newown == 1 ) then
            fedgex(nbndx) = clsx(i,j,k)
            iatm = atmx(i,j,k)
            aa = clsx(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm x3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - xi + REAL(i+1)
         xg(1) = gox + h*i + h*(1-aa)
         xg(2) = goy + h*j
         xg(3) = goz + h*k
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndx = nbndx + 1
            fedgex(nbndx) = 1.0d0 - aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
            iepsavx(4,nbndx) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
         nbndx = nbndx + 1
         aa = REAL(i+1) - xi - range3
         fedgex(nbndx) = 1.0d0 - aa
         if ( newown == 1 ) then
            fedgex(nbndx) = 1.d0 - clsx(i,j,k)
            iatm = atmx(i,j,k)
            aa = clsx(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm x4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
         iepsavx(4,nbndx) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = xi + range3 - REAL(i)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i + h*aa
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) =  aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(yi-j)**2
         if ( front >= 0.0d0 ) then
            range3 = sqrt(front)
            aa = xi - range3 - REAL(i)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) = aa
               if ( newown == 1 ) then
                  fedgex(nbndx) = clsx(i,j,k)
                  iatm = atmx(i,j,k)
                  aa = clsx(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm x5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(i+1) - xi + range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i + h*(1-aa)
            xg(2) = goy + h*j
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndx = nbndx + 1
               fedgex(nbndx) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(yi-j)**2
         if ( front >= 0.0d0 ) then
            range3 = sqrt(front)
            aa = REAL(i+1) - xi - range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndx = nbndx + 1
               fedgex(nbndx) = 1.0d0 - aa
               if ( newown == 1 ) then
                  fedgex(nbndx) = 1.d0 - clsx(i,j,k)
                  iatm = atmx(i,j,k)
                  aa = clsx(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm x6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
               iepsavx(4,nbndx) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(yi-j)**2
      nbndx = nbndx + 1
      range3 = sqrt(front)
      aa = REAL(i+1) - xi - range3
      fedgex(nbndx) = 1.0d0 - aa

      if ( newown == 1 ) then
         fedgex(nbndx) = 1.d0 - clsx(i,j,k)
         iatm = atmx(i,j,k)
         aa = clsx(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm x1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
      iepsavx(4,nbndx) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(yi-j)**2
      nbndx = nbndx + 1
      range3 = sqrt(front)
      aa = xi - range3 - REAL(i)
      fedgex(nbndx) = aa

      if ( newown == 1 ) then
         fedgex(nbndx) = clsx(i,j,k)
         iatm = atmx(i,j,k)
         aa = clsx(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm x2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavx(1,nbndx) = i; iepsavx(2,nbndx) = j; iepsavx(3,nbndx) = k
      iepsavx(4,nbndx) = -iatm
   end if

end subroutine epsfracx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = range3 + yi - REAL(j)
         fedgey(nbndy) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         nbndy = nbndy + 1
         aa = range3 - yi + REAL(j+1)
         fedgey(nbndy) = 1.0d0 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + yi - REAL(j)
         xg(1) = gox + h*i
         xg(2) = goy + h*j + h*aa
         xg(3) = goz + h*k
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndy = nbndy + 1
            fedgey(nbndy) = aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
            iepsavy(4,nbndy) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         nbndy = nbndy + 1
         aa = yi - range3 - REAL(j)
         fedgey(nbndy) = aa
         if ( newown == 1 ) then
            fedgey(nbndy) = clsy(i,j,k)
            iatm = atmy(i,j,k)
            aa = clsy(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm y3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - yi + REAL(j+1)
         xg(1) = gox + h*i
         xg(2) = goy + h*j + h*(1-aa)
         xg(3) = goz + h*k
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndy = nbndy + 1
            fedgey(nbndy) = 1.0d0 - aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
            iepsavy(4,nbndy) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
         nbndy = nbndy + 1
         aa = REAL(j+1) - yi - range3
         fedgey(nbndy) = 1.0d0 - aa
         if ( newown == 1 ) then
            fedgey(nbndy) = 1.d0 - clsy(i,j,k)
            iatm = atmy(i,j,k)
            aa = clsy(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm y4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
         iepsavy(4,nbndy) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(xi-i)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = yi + range3 - REAL(j)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j + h*aa
            xg(3) = goz + h*k
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(xi-i)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = yi - range3 - REAL(j)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = aa
               if ( newown == 1 ) then
                  fedgey(nbndy) = clsy(i,j,k)
                  iatm = atmy(i,j,k)
                  aa = clsy(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm y5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if 
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(zi-k)**2-(xi-i)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(j+1) - yi + range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j + h*(1-aa)
            xg(3) = goz + h*k
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndy = nbndy + 1
               fedgey(nbndy) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(zi-k)**2-(xi-i)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = REAL(j+1) - yi - range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndy = nbndy + 1
               fedgey(nbndy) = 1.0d0 - aa
               if ( newown == 1 ) then
                  fedgey(nbndy) = 1.d0 - clsy(i,j,k)
                  iatm = atmy(i,j,k)
                  aa = clsy(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm y6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
               iepsavy(4,nbndy) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(xi-i)**2
      nbndy = nbndy + 1
      range3 = sqrt(front)
      aa = REAL(j+1) - yi - range3
      fedgey(nbndy) = 1.0d0 - aa

      if ( newown == 1 ) then
         fedgey(nbndy) = 1.d0 - clsy(i,j,k)
         iatm = atmy(i,j,k)
         aa = clsy(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm y1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
      iepsavy(4,nbndy) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(zi-k)**2-(xi-i)**2
      nbndy = nbndy + 1
      range3 = sqrt(front)
      aa = yi - range3 - REAL(j)
      fedgey(nbndy) = aa

      if ( newown == 1 ) then
         fedgey(nbndy) = clsy(i,j,k)
         iatm = atmy(i,j,k)
         aa = clsy(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm y2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavy(1,nbndy) = i; iepsavy(2,nbndy) = j; iepsavy(3,nbndy) = k
      iepsavy(4,nbndy) = -iatm
   end if

end subroutine epsfracy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   integer flag, add_flag, flag_sub
   _REAL_ range1, range3, xi, yi, zi, aa, da, db
   _REAL_ front
   _REAL_ xg(3)

   if ( a == 2 .and. b == -2 ) then
      iatm = a1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = range3 + zi - REAL(k)
         fedgez(nbndz) = aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = iatm
      else
         epsint = depsout
      end if
   else if ( a == -2 .and. b == 2 ) then
      iatm = b1
      if ( sasopt > 0 ) then
         range1 = (radi(iatm)+dprob)*rh
      else if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         nbndz = nbndz + 1
         aa = range3 - zi + REAL(k+1)
         fedgez(nbndz) = 1.0d0 - aa
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = iatm
      else
         epsint = depsout
      end if
   else if ( a == 2 .and. b == -1 ) then
      flag_sub = 1
      iatm = a1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 + zi - REAL(k)
         xg(1) = gox + h*i
         xg(2) = goy + h*j
         xg(3) = goz + h*k + h*aa
         call flag_value(xg,b1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndz = nbndz + 1
            fedgez(nbndz) = aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
            iepsavz(4,nbndz) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = b1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
         nbndz = nbndz + 1
         aa = zi - range3 - REAL(k)
         fedgez(nbndz) = aa
         if ( newown == 1 ) then
            fedgez(nbndz) = clsz(i,j,k)
            iatm = atmz(i,j,k)
            aa = clsz(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm z3", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      end if
   else if ( b == 2 .and. a == -1 ) then
      flag_sub = 1
      iatm = b1
      if ( level == nfocus ) then
         range1 = radi(iatm)*rh
      else
         range1 = radip3(iatm)*rh
      end if
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
      if ( range3 > 0.0d0 ) then
         aa = range3 - zi + REAL(k+1)
         xg(1) = gox + h*i
         xg(2) = goy + h*j
         xg(3) = goz + h*k + h*(1-aa)
         call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
         if ( flag == 0 ) then ! the fraction does not lies in the reentry
            flag_sub = 0
            nbndz = nbndz + 1
            fedgez(nbndz) = 1.0d0 - aa
            epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
            iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
            iepsavz(4,nbndz) = iatm
         end if
      else
         flag_sub = 0
         epsint = epsout
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         range3 = sqrt(range1**2-(xi-i)**2-(yi-j)**2)
         nbndz = nbndz + 1
         aa = REAL(k+1) - zi - range3
         fedgez(nbndz) = 1.0d0 - aa
         if ( newown == 1 ) then
            fedgez(nbndz) = 1.d0 - clsz(i,j,k)
            iatm = atmz(i,j,k)
            aa = clsz(i,j,k)
            if ( iatm == 0 ) then
               write(6,'(a,3i6)') "wrong iatm z4", i, j, k
               stop
            end if
         end if
         epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
         iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
         iepsavz(4,nbndz) = -iatm
      end if
   else if ( a == 1 .and. b == -2 ) then
      flag_sub = 1
      iatm = b1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(xi-i)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = zi + range3 - REAL(k)
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j
            xg(3) = goz + h*k + h*aa
            call flag_value(xg,a1,flag) ! check whether the fraction lies in the reentry region
            if ( flag == 0 ) then ! the fraction does not lies in the reentry
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         iatm = a1
         range1 = dprob*rh
         xi = (arccrd(1,iatm) - gox)*rh
         yi = (arccrd(2,iatm) - goy)*rh
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(xi-i)**2-(yi-j)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = zi - range3 - REAL(k)
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = aa
               if ( newown == 1 ) then
                  fedgez(nbndz) = clsz(i,j,k)
                  iatm = atmz(i,j,k)
                  aa = clsz(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm z5", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( b == 1 .and. a == -2  ) then
      flag_sub = 1
      iatm = a1
      range1 = radi(iatm)*rh
      xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
      front = range1**2-(xi-i)**2-(yi-j)**2
      if ( front > 0.d0 ) then
         range3 = sqrt(front)
         aa = REAL(k+1) - zi + range3
         if ( aa > 0.d0 .and. aa < 1.d0 ) then 
            xg(1) = gox + h*i
            xg(2) = goy + h*j
            xg(3) = goz + h*k + h*(1-aa)
            call flag_value(xg,b1,flag)
            if ( flag == 0 ) then
               flag_sub = 0
               nbndz = nbndz + 1
               fedgez(nbndz) = 1.0d0 - aa
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         range1 = dprob*rh
         iatm = b1
         xi = (arccrd(1,iatm) - gox)*rh 
         yi = (arccrd(2,iatm) - goy)*rh 
         zi = (arccrd(3,iatm) - goz)*rh
         front = range1**2-(xi-i)**2-(yi-j)**2
         if ( front >= 0.d0 ) then
            range3 = sqrt(front)
            aa = REAL(k+1) - zi - range3
            if ( aa > 0.d0 .and. aa < 1.d0 ) then 
               flag_sub = 0 
               nbndz = nbndz + 1
               fedgez(nbndz) = 1.0d0 - aa
               if ( newown == 1 ) then
                  fedgez(nbndz) = 1.d0 - clsz(i,j,k)
                  iatm = atmz(i,j,k)
                  aa = clsz(i,j,k)
                  if ( iatm == 0 ) then
                     write(6,'(a,3i6)') "wrong iatm z6", i, j, k
                     stop
                  end if
               end if
               epsint = (depsin*depsout)/(depsin*(1.0d0-aa) + depsout*aa)
               iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
               iepsavz(4,nbndz) = -iatm
            end if
         end if
      end if
      if ( flag_sub == 1 ) then
         epsint = depsin
      end if
   else if ( a == -1 .and. b == 1 ) then 
      iatm = a1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(xi-i)**2-(yi-j)**2
      nbndz = nbndz + 1
      range3 = sqrt(front)
      aa = REAL(k+1) - zi - range3
      fedgez(nbndz) = 1.0d0 - aa

      if ( newown == 1 ) then
         fedgez(nbndz) = 1.d0 - clsz(i,j,k)
         iatm = atmz(i,j,k)
         aa = clsz(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm z1", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
      iepsavz(4,nbndz) = -iatm
   else if ( b == -1 .and. a == 1 ) then 
      iatm = b1
      range1 = dprob*rh
      xi = (arccrd(1,iatm) - gox)*rh
      yi = (arccrd(2,iatm) - goy)*rh
      zi = (arccrd(3,iatm) - goz)*rh
      front = range1**2-(xi-i)**2-(yi-j)**2
      nbndz = nbndz + 1
      range3 = sqrt(front)
      aa = zi - range3 - REAL(k)
      fedgez(nbndz) = aa

      if ( newown == 1 ) then
         fedgez(nbndz) = clsz(i,j,k)
         iatm = atmz(i,j,k)
         aa = clsz(i,j,k)
         if ( iatm == 0 ) then
            write(6,'(a,3i6)') "wrong iatm z2", i, j, k
            stop
         end if
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
      iepsavz(1,nbndz) = i; iepsavz(2,nbndz) = j; iepsavz(3,nbndz) = k
      iepsavz(4,nbndz) = -iatm
   end if

end subroutine epsfracz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine flag_value(xg,iatm,flag)

   use solvent_accessibility, only: ntri, triarc
   implicit none

   ! passed variables

   _REAL_ xg(3)
   integer iatm, flag

   ! local variables

   integer ii, jj, iarc
   _REAL_ xii(3), xjj(3), xij(3), dx(3), rxij, rd
   _REAL_ cosgij, cosgji, cosaij, cosaji, arcpos(3)
     
   arcpos(1:3) = arccrd(1:3,iatm)

   if ( iatm > narcdot - ntri ) then
      do ip = 1, 3
         iarc = triarc(ip,iatm-narcdot+ntri)
         jj = arcatm(1,iarc)
         xjj(1) = acrd(1,jj)
         xjj(2) = acrd(2,jj)
         xjj(3) = acrd(3,jj)
         ii = arcatm(2,iarc)
         xii(1) = acrd(1,ii)
         xii(2) = acrd(2,ii)
         xii(3) = acrd(3,ii)
         cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
         rxij = savarc(3,iarc)
         xij = rxij*(xjj - xii)

         dx = xg - xii; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

         dx = xg - xjj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
         cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

         if ( cosgij > cosaij .and. cosgji > cosaji ) then
            flag = 1
            exit
         else
            flag = 0
         end if
      end do
   else
      iarc = dotarc(iatm)

      jj = arcatm(1,iarc)
      xjj(1) = acrd(1,jj)
      xjj(2) = acrd(2,jj)
      xjj(3) = acrd(3,jj)
      ii = arcatm(2,iarc)
      xii(1) = acrd(1,ii)
      xii(2) = acrd(2,ii)
      xii(3) = acrd(3,ii)
      cosaij = savarc(1,iarc); cosaji = savarc(2,iarc)
      rxij = savarc(3,iarc)
      xij = rxij*(xjj - xii)

      dx = xg - xii; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
      cosgij =  (xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

      dx = xg - xjj; rd = 1.0d0/sqrt(dx(1)**2+dx(2)**2+dx(3)**2); dx = rd*dx
      cosgji = -(xij(1)*dx(1)+xij(2)*dx(2)+xij(3)*dx(3))!-small

      if ( cosgij > cosaij .and. cosgji > cosaji ) then
         flag = 1
      else
         flag = 0
      end if
   end if

end subroutine flag_value
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for x-edges with the level set function
subroutine epsfracx_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0

   if ( a > 0 ) then
      x1 = dble(i-1)
      x2 = dble(i  )
      x3 = dble(i+1)
      f1 = u(i-1,j,k)
      f2 = u(i  ,j,k)
      f3 = u(i+1,j,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i+1,j,k) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(i  )
      x2 = dble(i+1)
      x3 = dble(i+2)
      f1 = u(i  ,j,k)
      f2 = u(i+1,j,k)
      f3 = u(i+2,j,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(i)
   else
      aa = dble(i+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if ( sasopt == 2 ) nbndx = nbndx + 1
   fedgex(nbndx) = t - dble(i)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)

end subroutine epsfracx_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for y-edges with the level set function
subroutine epsfracy_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(j-1)
      x2 = dble(j  )
      x3 = dble(j+1)
      f1 = u(i,j-1,k)
      f2 = u(i,j  ,k)
      f3 = u(i,j+1,k)
      call root(x2,x3,x1,f2,f3,f1,t)
   end if

   ! if (i,j+1,k) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(j  )
      x2 = dble(j+1)
      x3 = dble(j+2)
      f1 = u(i,j  ,k)
      f2 = u(i,j+1,k)
      f3 = u(i,j+2,k)
      call root(x1,x2,x3,f1,f2,f3,t)
   end if

   if ( a > 0 ) then
      aa = t - dble(j)
   else
      aa = dble(j+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if ( sasopt == 2 ) nbndy = nbndy + 1
   fedgey(nbndy) = t - dble(j)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)

end subroutine epsfracy_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ assign fractional eps values for z-edges with the level set function
subroutine epsfracz_r( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout, u )

   implicit none

   ! passed variables

   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout
   _REAL_ u(0:xm+1,0:ym+1,0:zm+1)

   ! local variables

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa
   _REAL_ x1,x2,x3,f1,f2,f3,t

   ! if (i,j,k) is inside, f2 < 0 and f3 > 0 

   if ( a > 0 ) then
      x1 = dble(k-1)
      x2 = dble(k)
      x3 = dble(k+1)
      f1 = u(i,j,k-1)
      f2 = u(i,j,k)
      f3 = u(i,j,k+1)
      call root(x2,x3,x1,f2,f3,f1,t)
if ( i == 103 .and. j == 90 .and. k == 77 ) write(200, *) i,j,k, a, b, t
if ( i == 103 .and. j == 90 .and. k == 77 ) write(200, *) i,j,k, x2,x3,x1, f2,f3,f1
   end if

   ! if (i,j,k+1) is inside, f2 < 0 and f1 > 0

   if ( b > 0 ) then
      x1 = dble(k)
      x2 = dble(k+1)
      x3 = dble(k+2)
      f1 = u(i,j,k)
      f2 = u(i,j,k+1)
      f3 = u(i,j,k+2)
      call root(x1,x2,x3,f1,f2,f3,t)
if ( i == 103 .and. j == 90 .and. k == 77 ) write(200, *) i,j,k, a, b, t
if ( i == 103 .and. j == 90 .and. k == 77 ) write(200, *) i,j,k, x1,x2,x3, f1,f2,f3
   end if

   if ( a > 0 ) then
      aa = t - dble(k)
   else
      aa = dble(k+1) - t
   end if
   if ( abs(aa) < 1.d-12 ) aa = 0.d0
   if ( sasopt == 2 ) nbndz = nbndz + 1
   fedgez(nbndz) =  t - dble(k)
   epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
if ( i == 103 .and. j == 90 .and. k == 77 ) write(200, *) i,j,k, t, aa

end subroutine epsfracz_r
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ root returns the approximated root between x0 and x1 if f0*f1 <=0 using
!+ quadratic interpolation
subroutine root(x0,x1,x2,f0,f1,f2,t0)

   implicit none

   ! passed variables

   _REAL_ x0,x1,x2,f0,f1,f2,t0

   ! local variables

   _REAL_ b,c,a0,b0,c0,t,r1,r2

   b = (f0-f1)/(x0-x1)
   c = f2 - f1 - b*(x2-x1)
   c = c/( (x2-x0)*(x2-x1))

   a0 = c
   b0 = b - c*(x0+x1)
   c0 = f1 -b*x1 + c*x0*x1

   if ( a0 == 0 ) then
      t0 = -c0/b0
      return
   end if

   t = b0*b0 - 4.0d0*a0*c0

   ! If t <=0, must be double root t is close to zero

   if ( t <= 0.0d0 ) then
      t0 = -b0/(2.0d0*a0)
      return
   end if

   t = sqrt(t)
   if ( b0 >= 0.0d0 ) then
      r1 = (-b0-t)/(2.0d0*a0)
   else
      r1 = (-b0+t)/(2.0d0*a0)
   end if

   r2 = -b0/a0-r1

   if ( x0 <= r1 + 1.0d-7 .and. r1 <= x1+1.0d-7 ) then
      t0 = r1
   else
      t0 = r2
   end if

   if ( x0 > t0 ) t0 = x0
   if ( x1 < t0 ) t0 = x1

end subroutine root

end subroutine pb_exmol
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver for molecular dielectric map assignment.
subroutine pb_exmol_cap( pbverbose,ifcap )

   use poisson_boltzmann    
   use solvent_accessibility

   implicit none

   ! Passed variables
 
   logical pbverbose
   integer ifcap
 
   ! Local variables

   logical ses 
   integer ip, iatm, nwarn, xmymzm_ext
   _REAL_ xi, yi, zi
   _REAL_ range1, rh

   ! local array setup
   xmymzm_ext = xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8
 
   epsx(1:xmymzm+ym*zm) = epsout; epsy(1:xmymzm+xm*zm) = epsout; epsz(1:xmymzm+xm*ym) = epsout

   ! mark volume within vdw srf as 2, outside -2, so there are only contact-type
   ! boundary grid points

   rh = 1/h
   insas(1:xmymzm_ext) = -2
   zv(1:xmymzm_ext) = 9999.0d0
   iatm = -1
   range1 = radi(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   iatm = 1
   call exvwsph_cap( 2, insas, atmsas, zv )

   ! finally, save boundary edges for db energy and forces105 
   call epsbnd_cap( atmsas, insas )

   ! use the insas grid to setup epsx, epsy and epsz maps

   call epsmap_cap( insas, atmsas, epsx, epsy, epsz )


contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mark grid points within atomic vdw surf
subroutine exvwsph_cap( dielsph,insph,inatm,dst )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Mark grid points within an atom (dielectric constant dielsph) of index iatm
   ! as dielpsh. Modified from UHBD (Comp. Phys. Comm. 91:57-95, 1995) routines
   ! excrx() and exsrfx() by Michael Gilson and Malcom Davis.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none
    
   ! Passed variables
    
   integer  dielsph
   integer  insph(0:xm+1,0:ym+1,0:zm+1)
   integer  inatm(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ dst(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer  i, j, k
   integer  lowi, lowj, lowk
   integer  highi, highj, highk
   _REAL_ range2, range3, d2
    
   if ( zi+range1<0 .or. zi-range1>zm+1 ) return
   lowk = max(0,ceiling(zi - range1)); highk = min(zm+1,floor(zi + range1))
   do k = lowk, highk

      range2 = sqrt(range1**2-(zi-dble(k))**2)
      if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle
      lowj = max(0,ceiling(yi - range2)); highj = min(ym+1,floor(yi + range2))
      do j = lowj, highj

         range3 = sqrt(range2**2-(yi-dble(j))**2)
         if ( range3 >= 0.0d0 ) then

            if ( xi+range3<0 .or. xi-range3>xm+1 ) cycle
            lowi = max(0,ceiling(xi - range3)); highi = min(xm+1,floor(xi + range3))
            do i = lowi, highi
               d2 = (i-xi)**2 + (j-yi)**2 + (k-zi)**2
               if ( insph(i,j,k) == dielsph ) then
                  if ( d2 < dst(i,j,k) ) then
                     inatm(i,j,k) = iatm; dst(i,j,k) = d2
                  end if
                  cycle
               end if
               insph(i,j,k) = dielsph;
               inatm(i,j,k) = iatm; dst(i,j,k) = d2
            end do  ! i = lowi, highi
             
         end if  ! ( range3 > 0.0d0 )
          
      end do  ! j = lowj, highj
       
   end do  ! k = lowk, highk
    
end subroutine exvwsph_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Save dielectric boundary grid points
subroutine epsbnd_cap ( atmsas,insas )
    
   implicit none
    
   integer atmsas(0:xm+1,0:ym+1,0:zm+1)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   logical boundary
   integer buffer, i, j, k, clstmp
    
   nwarn = 0
   nbnd = 0
   buffer = 1
   do k = buffer, zm+1-buffer; do j = buffer, ym+1-buffer; do i = buffer, xm+1-buffer
       
      ! set up condition for a boundary grid point
       
      boundary = .false.
      if ( (insas(i,j,k) ==  1 .or. insas(i,j,k) ==  2) .and.&
           (insas(i-1,j,k) == -1 .or. insas(i-1,j,k) == -2 .or. insas(i+1,j,k) == -1 .or.&
            insas(i+1,j,k) == -2 .or. insas(i,j-1,k) == -1 .or. insas(i,j-1,k) == -2 .or.&
            insas(i,j+1,k) == -1 .or. insas(i,j+1,k) == -2 .or. insas(i,j,k-1) == -1 .or.&
            insas(i,j,k-1) == -2 .or. insas(i,j,k+1) == -1 .or. insas(i,j,k+1) == -2) ) then 
            boundary = .true.
      else if ( (insas(i,j,k) == -1 .or. insas(i,j,k) == -2) .and.&
           (insas(i-1,j,k) ==  1 .or. insas(i-1,j,k) ==  2 .or. insas(i+1,j,k) ==  1 .or.&
            insas(i+1,j,k) ==  2 .or. insas(i,j-1,k) ==  1 .or. insas(i,j-1,k) ==  2 .or.&
            insas(i,j+1,k) ==  1 .or. insas(i,j+1,k) ==  2 .or. insas(i,j,k-1) ==  1 .or.&
            insas(i,j,k-1) ==  2 .or. insas(i,j,k+1) ==  1 .or. insas(i,j,k+1) ==  2) ) then
            boundary = .true.
      end if
      if ( .not. boundary ) cycle
 
      nbnd = nbnd + 1; iepsav(1,nbnd) = i; iepsav(2,nbnd) = j; iepsav(3,nbnd) = k
      if ( ifcap /= 0 ) then
         iepsav(4,nbnd) = -1
         cycle
      end if
 
      ! for a grid point in contact region +/- 2 or in a  solvent probe, simply use the atom/probe that
      ! marks it
 
      clstmp = 0
      if ( abs(insas(i,j,k)) == 2 ) then
         clstmp = atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
      elseif ( insas(i,j,k) == -1 ) then
         clstmp = -atmsas(i,j,k)
         if ( clstmp == 0 ) then
            if ( pbverbose .and. level == nfocus ) write(6, '(a,4i4)') &
            'PB Warning in epsbnd(): No neighbor found for exposed boundary grid', i, j, k, insas(i,j,k)
            nwarn = nwarn + 1
         end if
 
      ! for a buried reentry grid point, find the atom that marked its neighoring exposed reentry
      ! grid points. Note that this may not be possible when grid spacing is large
 
      else if ( insas(i,j,k) == 1 ) then
         clstmp = -fndcls_cap( i, j, k, insas, atmsas )
         if ( clstmp == 0 ) then
            nwarn = nwarn + 1
         end if
      end if
 
      iepsav(4,nbnd) = clstmp
   end do; end do; end do
   if ( nwarn > 0 ) then
      if ( pbverbose .and. level == nfocus ) write(6, '(a,i4)') &
      'PB Warning in epsbnd(): No neighbor found for boundary grids total:', nwarn
   end if
 
end subroutine epsbnd_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the closest reentry probe for a reentry boundary grid
function fndcls_cap( i,j,k,insas,atmsas )
    
   implicit none
    
   ! Passed variables
    
   integer fndcls_cap, i, j, k
   integer  insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
    
   ! Local variables
    
   integer iatm, l, lp, ip, jp, kp, iip(6), jjp(6), kkp(6), clsatm(6)
   _REAL_ xg, yg, zg
   _REAL_ dx, dy, dz, d, clsdst, clscrd(3,6)

   ! first stack these candidates into a 1-d list
    
   iip(1)=i-1; iip(2)=i+1; jjp(1:2)=j; kkp(1:2)=k
   iip(3:4)=i; jjp(3)=j-1; jjp(4)=j+1; kkp(3:4)=k
   iip(5:6)=i; jjp(5:6)=j; kkp(5)=k-1; kkp(6)=k+1
   lp = 0
   do l = 1, 6
      ip = iip(l); jp = jjp(l); kp = kkp(l)
      if ( atmsas(ip,jp,kp) == 0 .or. insas(ip,jp,kp) /= -1 ) cycle
      lp = lp + 1; iatm = atmsas(ip,jp,kp); clsatm(lp) = iatm
      clscrd(1,lp) = arccrd(1,iatm)
      clscrd(2,lp) = arccrd(2,iatm)
      clscrd(3,lp) = arccrd(3,iatm)
   end do
 
   ! now find the closest
 
   xg = gox + i*h; yg = goy + j*h; zg = goz + k*h
   clsdst = 999.d0
   fndcls_cap = 0
   do ip = 1, lp
      dx = clscrd(1,ip) - xg; dy = clscrd(2,ip) - yg; dz = clscrd(3,ip) - zg
      d = abs(sqrt(dx**2 + dy**2 + dz**2) - dprob)
      if ( d >= clsdst ) cycle
      clsdst = d
      fndcls_cap = clsatm(ip)
   end do
 
end function fndcls_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map insas into epsmap
subroutine epsmap_cap( insas,atmsas,epsx,epsy,epsz )

   implicit none
   integer insas(0:xm+1,0:ym+1,0:zm+1), atmsas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)

   integer i, j, k, a, b, c, d, a1, b1, c1, d1
   _REAL_ epsint

   epsint = 2.0d0*epsin*epsout/(epsin+epsout)
 
   do k = 0, zm; do j = 0, ym; do i = 0, xm
      a = insas(i,j,k)
      b = insas(i+1,j,k)
      a1 = atmsas(i,j,k)
      b1 = atmsas(i+1,j,k)
      if ( j == 0 .or. k == 0 ) then
         ! do nothing
      else if ( sign(a,b) == a ) then
         if ( a > 0 ) then
            epsx(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracx_cap(i,j,k,a,b,a1,b1,rh,epsint,epsin,epsout)
         epsx(i,j,k) = epsint
      end if
      c = insas(i,j+1,k)
      c1 = atmsas(i,j+1,k)
      if ( i == 0 .or. k == 0 ) then
         ! do nothing
      else if ( sign(a,c) == a ) then
         if ( a > 0 ) then
            epsy(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracy_cap(i,j,k,a,c,a1,c1,rh,epsint,epsin,epsout)
         epsy(i,j,k) = epsint
      end if
      d = insas(i,j,k+1)
      d1 = atmsas(i,j,k+1)
      if ( i == 0 .or. j == 0 ) then
         ! do nothing
      else if ( sign(a,d) == a ) then
         if ( a > 0 ) then
            epsz(i,j,k) = epsin
         end if
      else
         if ( smoothopt == 1 ) call epsfracz_cap(i,j,k,a,d,a1,d1,rh,epsint,epsin,epsout)
         epsz(i,j,k) = epsint
      end if
   end do; end do; end do
 
!   do k = 1, zm
!      write(20, *) 'plane', k
!   do j = 1, ym
!      write(20, '(100f6.1)') epsx(1:xm,j,k)/epsin
!   end do
!   end do
!   do k = 1, zm
!      write(21, *) 'plane', k
!   do i = 1, xm
!      write(21, '(100f6.1)') epsy(i,1:ym,k)/epsin
!   end do
!   end do
!   do j = 1, ym
!      write(22, *) 'plane', j
!   do i = 1, xm
!      write(22, '(100f6.1)') epsz(i,j,1:zm)/epsin
!   end do
!   end do

end subroutine epsmap_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for x-edges
subroutine epsfracx_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this x-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracx_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(zi-k)**2-(yi-j)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - xi + dble(i+1)
      else
         aa = range3 + xi - dble(i)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracx_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for y-edges
subroutine epsfracy_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this y-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracy_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   ! obtain the position and radius of the atom (probe) in grid unit

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(zi-k)**2-(xi-i)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - yi + dble(j+1)
      else
         aa = range3 + yi - dble(j)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracy_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assign fractional eps values for z-edges
subroutine epsfracz_cap( i, j, k, a, b, a1, b1, rh, epsint, depsin, depsout )

   implicit none
   integer  i, j, k, a, b, a1, b1
   _REAL_ rh, epsint
   _REAL_ depsin, depsout

   integer iatm
   _REAL_ range1, range3, xi, yi, zi, aa

   ! mjhsieh: warning eliminator
   iatm = -1
   ! locate the atom that is crossing this z-edge, note the order iatm is assigned

   if ( a == 2 ) then
      iatm = a1
   else if ( b == 2 ) then
      iatm = b1
   else
      write(6,'(a)') 'PB Bomb in epsfracz_cap(): iatm not initialized.'
      call mexit(6,1)
   end if

   range1 = radip3(iatm)*rh
   xi = gcrd(1,iatm); yi = gcrd(2,iatm); zi = gcrd(3,iatm)
   range3 = sqrt(range1**2-(yi-j)**2-(xi-i)**2)
   if ( range3 > 0.0d0 ) then
      if ( b == 2 ) then
         aa = range3 - zi + dble(k+1)
      else
         aa = range3 + zi - dble(k)
      end if
      epsint = (depsout*depsin)/(depsin*(1.0d0-aa) + depsout*aa)
   else
      epsint = depsout
   end if

   ! other situations will not be considered and use default epsint value

end subroutine epsfracz_cap

end subroutine pb_exmol_cap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ surface area calculation routine 1
subroutine calc_sa1(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt)

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   implicit none

#  include "pb_constants.h"
#  include "../include/md.h"

   ! Passed variables
   _REAL_ acrd(3,*)
   integer xm,ym,zm,xmymzm,nbndx,nbndy,nbndz
   !integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   integer iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer sasopt

   ! Local variables

   integer i, j, k, iatm, ip
   integer cnt_dot, rnt_dot, dim_dot, tri_dot
   _REAL_ x(3), crd(3)
   _REAL_ rn(1:3), rsphere, dr, r1, r2, r3, h2, hh
   _REAL_ ds1, total_s1, ds2, total_s2
   _REAL_ ds, total_s , cnt_s, rnt_s, dim_s, tri_s
   _REAL_ dss, tss, rx, ry, rz, ess, e0, e1

   ! mjhsieh: warning eliminator
   rsphere = -1d0
   h2 = h*h
   hh = HALF*h

   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO

   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
!     crd(1) = gox + h*i + fedgex(ip)*h; crd(2) = goy + h*j; crd(3) = goz + h*k
      crd(1) = gox + h*i + hh; crd(2) = goy + h*j; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in calc_sa1(): cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(1))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     rx = ONE/rn(1)
!     dss = atan(rx*(rn(2)+hh)*(rn(3)+hh)/sqrt(rn(1)**2+(rn(2)+hh)**2+(rn(3)+hh)**2)) &
!         - atan(rx*(rn(2)+hh)*(rn(3)-hh)/sqrt(rn(1)**2+(rn(2)+hh)**2+(rn(3)-hh)**2)) &
!         - atan(rx*(rn(2)-hh)*(rn(3)+hh)/sqrt(rn(1)**2+(rn(2)-hh)**2+(rn(3)+hh)**2)) &
!         + atan(rx*(rn(2)-hh)*(rn(3)-hh)/sqrt(rn(1)**2+(rn(2)-hh)**2+(rn(3)-hh)**2)) 
!     dss = dss*rx*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndx

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
!     crd(1) = gox + h*i; crd(2) = goy + h*j + fedgey(ip)*h; crd(3) = goz + h*k
      crd(1) = gox + h*i; crd(2) = goy + h*j + hh; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in calc_sa1(): cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(2))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3 
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(2)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     ry = ONE/rn(2)
!     dss = atan(ry*(rn(1)+hh)*(rn(3)+hh)/sqrt(rn(2)**2+(rn(1)+hh)**2+(rn(3)+hh)**2)) &
!         - atan(ry*(rn(1)+hh)*(rn(3)-hh)/sqrt(rn(2)**2+(rn(1)+hh)**2+(rn(3)-hh)**2)) &
!         - atan(ry*(rn(1)-hh)*(rn(3)+hh)/sqrt(rn(2)**2+(rn(1)-hh)**2+(rn(3)+hh)**2)) &
!         + atan(ry*(rn(1)-hh)*(rn(3)-hh)/sqrt(rn(2)**2+(rn(1)-hh)**2+(rn(3)-hh)**2)) 
!     dss = dss*ry*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndy

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
!     crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+fedgez(ip)*h
      crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+hh

      if ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in calc_sa1(): cannot find owner of boundary grid points' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
      end if

      dr = abs(rn(3))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3 
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(3)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      total_s = total_s + ds

!     rz = ONE/rn(3)
!     dss = atan(rz*(rn(1)+hh)*(rn(2)+hh)/sqrt(rn(3)**2+(rn(1)+hh)**2+(rn(2)+hh)**2)) &
!         - atan(rz*(rn(1)+hh)*(rn(2)-hh)/sqrt(rn(3)**2+(rn(1)+hh)**2+(rn(2)-hh)**2)) &
!         - atan(rz*(rn(1)-hh)*(rn(2)+hh)/sqrt(rn(3)**2+(rn(1)-hh)**2+(rn(2)+hh)**2)) &
!         + atan(rz*(rn(1)-hh)*(rn(2)-hh)/sqrt(rn(3)**2+(rn(1)-hh)**2+(rn(2)-hh)**2)) 
!     dss = dss*rz*rsphere**2*dr
!     tss = tss + dss
!     ess = abs(ds*h2-dss)/dss
!     e1 = ess + e1
!     if ( ess > e0 ) e0 = ess

      if ( iatm > 0 ) then
         cnt_s = cnt_s + ds
         cnt_dot = cnt_dot + 1
      else
         rnt_s = rnt_s + ds
         rnt_dot = rnt_dot + 1
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + ds
            tri_dot = tri_dot + 1
         else
            dim_s = dim_s + ds
            dim_dot = dim_dot + 1
         end if 
      end if

   end do !nbndz

   total_s = total_s*h2
!  total_s1 = total_s1*h2
!  total_s2 = total_s2*h2
   cnt_s = cnt_s*h2
   rnt_s = rnt_s*h2
   tri_s = tri_s*h2
   dim_s = dim_s*h2
!  e1 = e1/(nbndx+nbndy+nbndz)
!  if ( saopt < 0 .or. imin == 6 ) then
      write(6,'(a,f12.4)') ' Total molecular surface',total_s
!     write(6,'(a,f20.4)') 'Total contact surface',cnt_s
!     write(6,'(a,f20.4)') 'Total reentry surface',rnt_s
!     write(6,'(a,f20.4)') 'Total dimer surface',dim_s
!     write(6,'(a,f20.4)') 'Total trimer surface',tri_s
!     write(6,'(a,i10)') 'contact boundary point',cnt_dot
!     write(6,'(a,i10)') 'reentrant boundary point',rnt_dot
!     write(6,'(a,i10)') 'dimer boundary point',dim_dot
!     write(6,'(a,i10)') 'trimer boundary point',tri_dot

!     write(6,'(a,f20.4)') 'zero order',total_s1
!     write(6,'(a,f20.4)') 'second order',total_s2
!     write(6,'(a,f20.4)') 'all order',tss
!     write(6,'(a,f20.4)') 'mode oo error',e0
!     write(6,'(a,f20.4)') 'mode 1 error',e1
!     print *,total_s,cnt_s,rnt_s
!  end if

end subroutine calc_sa1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ surface area calculation routine 2
subroutine calc_sa2(acrd,xm,ym,zm,xmymzm,nbnd,iepsav,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,sasopt,smoothopt,&
                    epsin,epsout,insas,epsx,epsy,epsz)

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   implicit none

#  include "pb_constants.h"
#  include "../include/md.h"

   ! Passed variables

   _REAL_ acrd(3,*)
   integer xm,ym,zm,xmymzm,nbnd
   integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   integer sasopt,smoothopt
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ epsin, epsout
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)

   ! Local variables

   integer ip, i, j, k, iatm
   _REAL_ g(3), x(3), dx(3), hh(3), repsp(3), repsm(3)
   _REAL_ dist, dx2, ds, total_s, sgn, half_h, epsth, repsin, repsout
   _REAL_ cnt_s, rnt_s, term, dim_s, tri_s

   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   repsin = ONE/epsin; repsout = ONE/epsout
   epsth = TWO/(repsin+repsout)
!  half_h = HALF*h

   do ip = 1, nbnd

      ! collecting boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the surface grid point on to the molecular surface, crd() is the
      ! new coord, and x() is the atom/probe coord, fx/y/z0 is the grid version
      ! of crd()

      if      ( iatm == 0 ) then
         write(6,'(a)') 'PB Bomb in calc_sa2(): can not find owner of boundary grid points'
         call mexit(6, 1)
      else if ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            dist = radi(iatm)+dprob
         else
            dist = radi(iatm)
         end if
         sgn = ONE
      else
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         dist = dprob
         sgn = -ONE
      end if

      dx = g - x; dx2 = dx(1)**2 + dx(2)**2 + dx(3)**2
!     hh = half_h

      dist = sqrt(dx2)
!     ds = (dx(1)+hh(1))/epsx(i,j,k)-(dx(1)-hh(1))/epsx(i-1,j,k) + &
!          (dx(2)+hh(2))/epsy(i,j,k)-(dx(2)-hh(2))/epsy(i,j-1,k) + &
!          (dx(3)+hh(3))/epsz(i,j,k)-(dx(3)-hh(3))/epsz(i,j,k-1)  
!     ds = dx(1)/epsx(i,j,k)-dx(1)/epsx(i-1,j,k) + &
!          dx(2)/epsy(i,j,k)-dx(2)/epsy(i,j-1,k) + &
!          dx(3)/epsz(i,j,k)-dx(3)/epsz(i,j,k-1) 
      if ( smoothopt == 1 .or. smoothopt == 2 ) then
         repsp = repsin
         repsm = repsin
         if ( epsx(i,j,k) > epsth ) repsp(1) = repsout
         if ( epsy(i,j,k) > epsth ) repsp(2) = repsout
         if ( epsz(i,j,k) > epsth ) repsp(3) = repsout
         if ( epsx(i-1,j,k) > epsth ) repsm(1) = repsout
         if ( epsy(i,j-1,k) > epsth ) repsm(2) = repsout
         if ( epsz(i,j,k-1) > epsth ) repsm(3) = repsout
         ds = dx(1)*repsp(1)-dx(1)*repsm(1) + &
              dx(2)*repsp(2)-dx(2)*repsm(2) + &
              dx(3)*repsp(3)-dx(3)*repsm(3) 
      else
         ds = dx(1)/epsx(i,j,k)-dx(1)/epsx(i-1,j,k) + &
              dx(2)/epsy(i,j,k)-dx(2)/epsy(i,j-1,k) + &
              dx(3)/epsz(i,j,k)-dx(3)/epsz(i,j,k-1) 
      end if
      term = sgn*ds/dist
      total_s = total_s + term
      if ( iatm > 0  ) then
         cnt_s = cnt_s + term
      else
         rnt_s = rnt_s + term
         if ( -iatm > narcdot - ntri ) then
            tri_s = tri_s + term
         else
            dim_s = dim_s + term
         end if
      end if
      
   end do

   term = h*h/(repsout-repsin)
   total_s = total_s*term
   cnt_s = cnt_s*term
   rnt_s = rnt_s*term
   dim_s = dim_s*term
   tri_s = tri_s*term
!  if ( saopt < 0 .or. imin == 6 ) then
!     print *, nbnd, nbndx+nbndy+nbndz
      write(6,'(1x,a,f12.4)') 'Total molecular surface',total_s
!     write(6,'(1x,a,f12.4)') 'Total contact surface',cnt_s
!     write(6,'(1x,a,f12.4)') 'Total reentry surface',rnt_s
!     write(6,'(1x,a,f12.4)') 'Total dimer surface',dim_s
!     write(6,'(1x,a,f12.4)') 'Total trimer surface',tri_s
!     print *,total_s,cnt_s,rnt_s
!  end if

end subroutine calc_sa2
