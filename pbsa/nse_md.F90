#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

module nse_md

integer inse

_REAL_ rhow_effect,dprob
_REAL_,allocatable :: radi(:)

_REAL_ fillratio
integer xm,ym,zm
_REAL_ h,gox,goy,goz
_REAL_,allocatable :: lvlset(:,:,:)
 
! these should be allocatable
_REAL_ fx(1:xm,1:ym,1:zm), fy(1:xm,1:ym,1:zm), fz(1:xm,1:ym,1:zm)

contains

subroutine nse_read

   namelist /nse/ fillratio,dprob,h,rhow_effect
   integer ifind

   fillratio=2.0
   dprob=1.4
   h=0.5
   rhow_effect=1.129d0

   call nmlsrc('nse',5,ifind)
   if ( ifind/=0 ) then
      read(5,nml=nse)
   end if

end subroutine nse_read
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ initialize the level set function using the revised density approach
subroutine nse_init(verbose,pbprint,natom,ntypes,iac,ico,cn1,cn2,acrd)
    
   implicit none
    
#  include "pb_constants.h"
#  include "extra.h"
   ! Passed variables

   logical verbose, pbprint
   integer natom, nres, ntypes, iac(*), ico(*)
   _REAL_ cn1(*), cn2(*)
   _REAL_ acrd(3,natom)

   ! Local variables
    
   integer iatm, ic, ier
   _REAL_ xlength, ylength, zlength
   _REAL_ xbox, ybox, zbox
   _REAL_ xmin,xmax,ymin,ymax,zmin,zmax
   _REAL_ gcrd(3,natom)
   _REAL_ tmp

   ! part 1. Set Atomic Radii as Rmin

   if ( allocated(radi) ) then
      deallocate(radi, stat = ier); REQUIRE(ier==0)
   end if
   allocate(radi(1:natom), stat = ier); REQUIRE(ier==0)

   do iatm = 1, natom
      ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
      if (cn2(ic) /= ZERO) then
         tmp = (cn1(ic)/cn2(ic))**(SIXTH)/2 ! this is sigma
         radi(iatm) = tmp*(TWO**(SIXTH)) ! this is Rmin 
      else
         radi(iatm) = ZERO
      endif
   end do

   ! part 2. Set Finite Difference Grid
   ! first find bounding box
    
   if ( verbose .and. pbprint ) then
      write(6,'(a)')
      write(6,'(a)')
      write(6,'(a)') '======== Setting up NSE Grid ========'
      write(6,'(a)') 'Using bounding box for grid setup'
   end if

   xmin =  9999.0; ymin =  9999.0; zmin =  9999.0
   xmax = -9999.0; ymax = -9999.0; zmax = -9999.0
   do iatm = 1, natom 
      if ( acrd(1,iatm)-radi(iatm) .lt. xmin ) xmin = acrd(1,iatm)-radi(iatm)
      if ( acrd(1,iatm)+radi(iatm) .gt. xmax ) xmax = acrd(1,iatm)+radi(iatm)
      if ( acrd(2,iatm)-radi(iatm) .lt. ymin ) ymin = acrd(2,iatm)-radi(iatm)
      if ( acrd(2,iatm)+radi(iatm) .gt. ymax ) ymax = acrd(2,iatm)+radi(iatm)
      if ( acrd(3,iatm)-radi(iatm) .lt. zmin ) zmin = acrd(3,iatm)-radi(iatm)
      if ( acrd(3,iatm)+radi(iatm) .gt. zmax ) zmax = acrd(3,iatm)+radi(iatm)
   end do

   ! next find the box center, round it to the nearest h unit for easy restarting
   ! and the box dimension

   xbox = (xmax + xmin)/TWO; ybox = (ymax + ymin)/TWO; zbox = (zmax + zmin)/TWO
   xbox = nint(xbox/h)*h; ybox = nint(ybox/h)*h; zbox = nint(zbox/h)*h
   xlength = xmax-xmin; ylength = ymax-ymin; zlength = zmax-zmin
   if ( verbose .and. pbprint ) then
      write(6, '(1x,a,3f10.3)') 'Bounding Box Center:  ', xbox, ybox, zbox
      write(6, '(1x,a,3f10.3)') 'Xmin, Xmax, Xmax-Xmin:', xmin, xmax, xmax-xmin
      write(6, '(1x,a,3f10.3)') 'Ymin, Ymax, Ymax-Ymin:', ymin, ymax, ymax-ymin
      write(6, '(1x,a,3f10.3)') 'Zmin, Zmax, Zmax-Zmin:', zmin, zmax, zmax-zmin
   end if

   ! finally set the grid dimensions and the grid origin

   xm = nint( xlength*fillratio/h )
   ym = nint( ylength*fillratio/h )
   zm = nint( zlength*fillratio/h )
   gox = - dble(xm+1) * h * HALF + xbox
   goy = - dble(ym+1) * h * HALF + ybox
   goz = - dble(zm+1) * h * HALF + zbox
   if ( verbose .and. pbprint ) then
      write(6, '(a,i5,1x,3i5)') ' Grid dimension  ',  xm, ym, zm
      write(6, '(a,1x,3f10.3)') ' Grid origin  ', gox, goy, goz
      write(6,*) fillratio,h,xm,ym,zm,gox,goy,goz,natom,gcrd,radi(1)+2.0*dprob
   end if

   ! part 3. Convert Lab Frame Coords to Grid Frame 

   do iatm = 1, natom
      gcrd(1,iatm) = (acrd(1,iatm) - gox)/h
      gcrd(2,iatm) = (acrd(2,iatm) - goy)/h
      gcrd(3,iatm) = (acrd(3,iatm) - goz)/h
   end do

   ! part 4. Set up Level Set and Do Some Initialization

   call density_lvlset(xm,ym,zm,dprob*h,natom,gcrd)
   !call nse_init_calc() ! turn it on when it's ready


end subroutine nse_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ set level set function using the revised density approach
subroutine density_lvlset(xm,ym,zm,dprobh,natom,gcrd)

   ! WMSB: Add atomic density contribution using the inkblot method as in SES/SAS

   implicit none
   integer natom
   _REAL_ dprobh
   _REAL_ gcrd(3,natom)
   _REAL_ lvlset(0:xm+1,0:ym+1,0:zm+1)

   integer i, j, k, l, ier
   integer lowi, lowj, lowk
   integer highi, highj, highk
   integer xm,ym,zm
   _REAL_ xi, yi, zi, rh
   _REAL_ range0, range1, range2, range3
   _REAL_ dist

   rh = 1.0d0/h
   
   do l = 1, natom
      range0 = radi(l)*rh ! atom radius in gridpoint units
      if ( range0 == 0.0d0 ) cycle ! don't waste time on zero radius atoms
      xi = gcrd(1,l); yi = gcrd(2,l); zi = gcrd(3,l)
      range1 = range0 + dprobh * 2.d0 ! distance till atomic lvlset contrib -> 0 
      if ( zi+range1<0 .or. zi-range1>zm+1 ) cycle ! this shouldn't happen...
      lowk = max(1,ceiling(zi - range1)); highk = min(zm,floor(zi + range1))
      do k = lowk, highk ! z indices (grid line)    
         range2 = sqrt(range1**2-(zi-dble(k))**2)
         if ( yi+range2<0 .or. yi-range2>ym+1 ) cycle ! this shouldn't happen...
         lowj = max(1,ceiling(yi - range2)); highj = min(ym,floor(yi + range2))
         do j = lowj, highj ! y indices (grid disc)
            range3 = sqrt(range2**2-(yi-dble(j))**2)
            if ( range3 > 0.0d0 ) then ! sanity check on range3
               lowi = max(1,ceiling(xi-range3)); highi = min(xm,floor(xi+range3))
               do i = lowi, highi ! x indices (grid sphere)
                  dist = sqrt((xi-i)**2+(yi-j)**2+(zi-k)**2) - range0
                  dist = dist / (2.d0*dprobh)
                  lvlset(i,j,k) = lvlset(i,j,k) + density_atom(dist)
               end do ! loop over x indices (grid sphere)
            end if ! sanity check on range3
         end do ! loop over y indicies (grid disc)
      end do ! loop over z indicies (grid line)
   end do
   lvlset=1.0d0-lvlset

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Single atom revised density function
_REAL_ function density_atom(dist)

   implicit none

   _REAL_ dist

   integer m
   _REAL_, parameter :: dash(6) = (/0.00d0,0.20d0,0.40d0,0.60d0,0.80d0,1.00d0 /)
   _REAL_, parameter :: spcoef(5,4) = &
       reshape((/1.000000  ,0.2100000  ,0.1500000  ,0.0500000  ,0.010000   ,  &
                -4.527143  ,-2.067608  ,0.0475730  ,-0.522686  ,-0.056828  ,  &
                -3.640532  ,15.938209  ,-5.362303  ,2.5110050  ,-0.181716  ,  &
                32.631235  ,-35.500854 ,13.122180  ,-4.487867  ,1.079289/), (/5,4/))

   ! The distance is in the unit of solvent probe diameter

   density_atom = 0.0d0
   if ( dist > 1.d0 ) then
   else if ( dist < 0.d0 ) then
      density_atom = 1.0d0 - 4.527143d0 * dist
   else
      do m = 1, 5
         if ( dist > dash(m) .and. dist <= dash(m+1) ) then
            density_atom = density_atom   + &
            spcoef(m,1)                   + &
            spcoef(m,2)*(dist-dash(m))    + &
            spcoef(m,3)*(dist-dash(m))**2 + &
            spcoef(m,4)*(dist-dash(m))**3
         end if
      end do
   end if

end function density_atom


end subroutine density_lvlset
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculate van der Waals forces between atoms and fluid
subroutine nse_npforce(natom,ntypes,iac,ico,cn1,cn2,xatm,f,nstep)

   implicit none
   ! should also output energy in addition to force.
   ! turn on when ready
   !use subroutines, only: t,fx,fy,fz,heavyside,dx
   _REAL_ t, dx
  
#  include "extra.h"
#  include "pb_constants.h"

   integer natom, nres, ntypes,  iac(*), ico(*)
   _REAL_ cn1(*), cn2(*), xatm(3,natom), f(3,natom)
   _REAL_ enbrfcav, enbrfdis,cutoff,dist2,range1,range2,range3,heavy
   integer lowk,highk,lowj,highj,lowi,highi,nstep
    
   ! Local variables
   integer i,j,k
   integer ic, iatm

   _REAL_ mdsig_ow, mdsig_iatm, epsln_iatm, sigow2, sigow4, sigow6,cn1ij,cn2ij
   _REAL_ a(natom), b(natom), epsow(natom), sigow(natom), rminow(natom)
   _REAL_, parameter :: rhow = 3.3330d-02
  
   _REAL_ dxx,dyy,dzz,xii,yii,zii,d2inv,f2,f1,df,r6,fw1,fw2,fw3,ene,ene1

   mdsig_ow = 1.7683d0*(TWO**(-SIXTH))
   f=0.0
   ene=0.0
   ene1=0.0
   cutoff=6.0

   do iatm = 1, natom
      ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
      if (cn2(ic) /= ZERO) then
         mdsig_iatm = (cn1(ic)/cn2(ic))**SIXTH/2
         epsln_iatm = cn2(ic)/(256.0d0*mdsig_iatm**6)   
      else
         mdsig_iatm = ZERO
         epsln_iatm = ZERO
      endif
      sigow(iatm) = mdsig_iatm + mdsig_ow
      rminow(iatm) = sigow(iatm)*(TWO**SIXTH)
      epsow(iatm) = sqrt(epsln_iatm*0.1520d0)
      sigow2 = sigow(iatm)*sigow(iatm)
      sigow4 = sigow2*sigow2
      sigow6 = sigow2*sigow4
      b(iatm) = FOUR*epsow(iatm)*sigow6*rhow_effect
      a(iatm) = b(iatm)*sigow6
      xii = (xatm(1,iatm)-gox)/h; yii = (xatm(2,iatm)-goy)/h; zii =(xatm(3,iatm)-goz)/h
      range1=cutoff/h

      lowk = ceiling(zii - range1); highk = floor(zii + range1)
      do k = lowk, highk ! z indices (grid line)    
         range2 = sqrt(range1**2-(zii-dble(k))**2)
         lowj = ceiling(yii - range2); highj = floor(yii + range2)
         do j = lowj, highj ! y indices (grid disc)
            range3 = sqrt(range2**2-(yii-dble(j))**2)
            lowi = ceiling(xii-range3); highi = floor(xii+range3)
            do i = lowi, highi ! x indices (grid sphere)
               if (i<=1 .or. j<=1 .or.k<=1 .or. i>=xm-1 .or. j>=ym-1 .or. k>=zm-1)then
                  heavy=1.0d0
               else
                  heavy=1.0d0!-heavyside(t(i,j,k),dx) ! turn on when ready
               end if
               if (heavy==0.0) cycle
               dxx = xatm(1,iatm) - (gox+i*h); dyy = xatm(2,iatm) - (goy+j*h); dzz = xatm(3,iatm) -(goz+k*h) ;
                                                
               d2inv = ONE/(dxx**2+dyy**2+dzz**2 )
               cn1ij = a(iatm) ; cn2ij = b(iatm)
               r6 = d2inv**3; f2 = cn2ij*r6; f1 = cn1ij*(r6*r6)
               df = ( SIX*( (f2-f1)-f1 ) )*d2inv
               ene=ene+(f1-f2)*(h**3)*rhow*heavy
               df = df*rhow*heavy
                      
               fw1 = dxx*df; fw2 = dyy*df; fw3 = dzz*df
               f(1,iatm) = f(1,iatm) - fw1*(h**3)
               f(2,iatm) = f(2,iatm) - fw2*(h**3)
               f(3,iatm) = f(3,iatm) - fw3*(h**3)
          
               if(i<=1 .or. j<=1 .or.k<=1 .or. i>=xm-1 .or. j>=ym-1 .or. k>=zm-1) cycle
          
               fx(i,j,k) = fx(i,j,k) + fw1*6.95d19
               fy(i,j,k) = fy(i,j,k) + fw2*6.95d19
               fz(i,j,k) = fz(i,j,k) + fw3*6.95d19
            end do
        end do
     end do
     ene=ene-4*pi*rhow*(-a(iatm)/9.0/(cutoff**9)+b(iatm)/3.0/(cutoff**3))
     ene1=ene1-4*pi*rhow*(-a(iatm)/9.0/(cutoff**9)+b(iatm)/3.0/(cutoff**3))
   end do

   !write(6,*) h,ene,ene1
   !rite(12,*) -4*pi*rhow*(-a(1)/9.0/(sigow(1)**9)+b(1)/3.0/(sigow(1)**3))
   do iatm = 1, natom
      !write(6,*) sigow(iatm),xatm(1,iatm),xatm(2,iatm),xatm(3,iatm)
      write(6,*) nstep,f(1,iatm),f(2,iatm),f(3,iatm),a(iatm),b(iatm)
   end do


end subroutine nse_npforce
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ do one step of time integration of nse
subroutine nse_onestep(temp0, nstep)

   _REAL_ temp0
   integer nstep


end subroutine nse_onestep


end module nse_md
