#include "../include/dprec.fh"
#include "copyright.h"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of module navier_stokes here]
module nse_time
   _REAL_ nsetime_bulletin (20)
   logical nsetimerswitch(20)
   character(len=25) NSETIME_desc (20)

contains


subroutine nsetimer_init
   implicit none
   nsetimerswitch=.false.
   nsetime_bulletin=0d0
   nsetime_desc(1) = 'regular update right-hand side'
   nsetime_desc(2)= 'IIM'
   nsetime_desc(3)='iccg'
   nsetime_desc(4)='compute interface pressure'
   nsetime_desc(5)='interpolate boundary condition'
   nsetime_desc(6)='compute residue'
   nsetime_desc(7)='total time'
   nsetime_desc(8)='total time for compute matrix'
   nsetime_desc(9)='total time for solve matrix'  
   nsetime_desc(10)='total time for non-parallel'  
end subroutine nsetimer_init

subroutine nsetimer_start( num_event )
   implicit none

   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if (nsetime_bulletin(num_event) <= 0) then
      nsetime_bulletin(num_event) = mytime
   else
      nsetime_bulletin(num_event) = mytime-nsetime_bulletin(num_event)
   endif
   nsetimerswitch(num_event) = .true.
   return
end subroutine nsetimer_start

subroutine nsetimer_stop( num_event )
   implicit none
   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if ( .not. nsetimerswitch(num_event) ) then
      write(6,*)'timer: need to start the time in nse before stop',num_event
      call mexit(6,1)
   end if
   nsetime_bulletin(num_event) = mytime-nsetime_bulletin(num_event)
   nsetimerswitch(num_event) = .false.
   return
end subroutine nsetimer_stop


subroutine nsetimer_summary(num_time_event)
   implicit none
   integer i,num_time_event

   do i = 1, num_time_event
      if ( .not. nsetimerswitch(i) ) then
         if (nsetime_bulletin(i) > 0) &
            write(6,"('|',1x,a,f10.2)") nsetime_desc(i),nsetime_bulletin(i)
      else
         write(6,*) "Warning, this timer is not stopped: ",nsetime_desc(i)
      endif
   enddo
   return
end subroutine nsetimer_summary
end module nse_time

module navier_stokes
   use nse_time
   use linp
   use aug_solver
   implicit none
   ! Following are used in parrallel
   !       use fishpack
 !  integer numtasks,mytaskid
#include "extra.h"
#include "files.h"
#ifdef MPI
   include "mpif.h"
#  include "parallel.h"
#  undef  MPI_MAX_PROCESSORS
#  define MPI_MAX_PROCESSORS 256
#endif
   !Common block is for compatability:
!   common/parallel/numtasks,mytaskid,notdone, &
!         iparpt,iparpt3,&
!         rcvcnt,rcvcnt3,&
!         mpi_orig
!    common/extra_logical/master

   ! nse_debug
   ! 0: normal run
   ! 1: checking interpx using analytical pressure for single ion
   ! 2. checking IIM using analytical jump conditions and sources for single ion
   integer npfjump,npinterp
   parameter(npfjump=7,npinterp=11) 
   integer nse_debug,velocity_solver,solverorder,numbubble
   logical ifupdate,ifsetpressure,ifzero
    _REAL_ solvernorm,nsernu,ratio_t,ctfactor(1:3),zonexyz(1:3,1:3)
    _REAL_,allocatable :: fjumpinfo(:,:)
    _REAL_,allocatable :: interpinfo(:,:)
    integer,allocatable :: numfjump(:),numinterp(:)
   ! parameters for the analytical test case

   _REAL_ x0
   _REAL_ y0
   _REAL_ z0
   _REAL_ r0

   contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nse here]
   subroutine nse(m, n, l, h, gox, goy, goz)
      implicit none
      ! passed variables from Amber
      integer  m,n,l
      _REAL_  gox,goy,goz
      _REAL_  h,xa,xb,ya,yb,za,zb
      !       parameter(m=128,n=128,l=128,n1max=10000,xa=-2.0d0,xb=2.0d0,
      !       parameter(m=96,n=96,l=96,n1max=10000,xa=-3.0d0,xb=3.0d0,
      !parameter(m=48,n=48,l=48,n1max=10000,xa=-2.0d0,xb=2.0d0,&
      !          ya=-2.0d0,yb=2.0d0,za=-2.0d0,zb=2.0d0,np=20,mm=50)
      integer  np
      parameter(np=22)
      !    1           ya=-3.0d0,yb=3.0d0,za=-3.0d0,zb=3.0d0,np=20,mm=50)

      !module variables
      !x,y,z: Coordinate
      !u,v,w: Velocity
      !p1,p: pressure
      !phi: level set
      !aret: domain volume
      !indexx,indexy,indexz: irregular point on (x,y,z)
      !cinfox,cinfoy,cinfoz: coordinate, normal and tangetial
      !direction of irregular point
      !index2: projection point
      !cinfo: coordinate,normal and tangetial direction of irregular
      !points
      !ft,fn,fm: jump condition on un,vn,wn
      !ftn: merge ft,fn,fm
      !bf: right-hand side
      !fvec: matrix coordinate
      !un: velocity on normal direction
      !kzone: domain on grid points
      !gmat:coefficient matrix of augmented variable

      _REAL_  x(0:m),y(0:n),z(0:l) 
      _REAL_  u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) 
      _REAL_  fw(0:m,0:n,0:l),gw(0:m,0:n,0:l),hw(0:m,0:n,0:l) 
      _REAL_  ui(0:m,0:n,0:l),vi(0:m,0:n,0:l),wi(0:m,0:n,0:l) 
      _REAL_  u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) 
      _REAL_  p1(0:m,0:n,0:l),p(0:m,0:n,0:l) 
      _REAL_  phi(0:m,0:n,0:l),phi0(0:m,0:n,0:l) 
      _REAL_  aret(100) 
      integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ un(0:m,0:n,0:l) 
      integer isign(0:m,0:n,0:l)
      _REAL_  snorm(m,n,l,2) 
      _REAL_ haja(0:m,0:n,0:l) 
      _REAL_ work(2,0:m,0:n,0:l) 
      integer kzone(0:m,0:n,0:l) 
      _REAL_ al,coc1,coc2,coc3,dle,dt, e1,e2,e3,e4,eng,et,hx,hy,hz
      integer i,j,k,i1,j1,k1,ki,kll,kmax,kout,kout1,n2,n22,ierr,ii
      integer nb,num,n3,nzone,mm,izonex,izoney,izonez,nlp,nep,kdm
       _REAL_ phx,phy,phz,pt,pt_std,restart,rnu,width,width2,sc,t,tfinal,rec_value
       _REAL_ tol_pt,tol_pt2,umax,unmax,unmin,unorm,tolphi,dt0,t_standard 

      _REAL_,allocatable :: cinfox(:),cinfoy(:),cinfoz(:),cinfo(:)&
                            ,fnt(:),bf(:),fvec(:)    

      real*4 tmp(2),time, dtime
      character(10) fname,fname1,fname2,fname3
      character(6) str
      ! setting parameters
      !rnu: viscosity
      !r0:radius of static boundary
      !hx,hy,hz: length of grid at the direction of x,y,z
      !h: minimim of hx,hy,hz
      !duration of time between each step
!      ! common /rnu/rnu
!      ! common /radius/r0
!      ! common /aret/aret
      !write(*,*) 'node',master,mytaskid
      solverorder=1
      numbubble=1
      t_standard=1.0d0/30.0d0
      kout = 1
      kout1 = 1
      nlp=10
      nep=4
      ! h is now passed in
      hx = h
      hy = h
      hz = h
      !hx = (xb-xa)/m
      !hy = (yb-ya)/n
      !hz = (zb-za)/l
      !h = min(hx,hy,hz)

      ! xa, xb etc are now computed from h, l, m, n, and grid origin
      xa = gox + 1*h
      ya = goy + 1*h
      za = goz + 1*h
      xb = gox + (m+1)*h
      yb = goy + (n+1)*h
      zb = goz + (l+1)*h
      ! r0 is now passed in
      !r0 = 0.8d0

      rnu = nsernu

      dt = t_standard
      ratio_t=1.0d0
      kmax =1
      t = 0.0d0
      tfinal = 100.1d0/3.0d0

      kmax=int(tfinal/dt)

      do i=0,m
         x(i) = xa + hx*float(i)
      end do

      do j=0,n
         y(j) = ya + hy*float(j)
      end do

      do k=0,l
         z(k) = za + hz*float(k)
      end do
      call nsetimer_init()         

      !initialize solution vectors
      !ue,ve,we,pe doesn't depend on t
      sc = 1.0d0
               if(nse_debug .le. 8 .and. nse_debug .ne.0) then
      do k=0,l
         do j=0,n
            do i=0,m
               u1(i,j,k)= ue(t-dt,rnu,x(i),y(j),z(k))/sc
               v1(i,j,k)= ve(t-dt,rnu,x(i),y(j),z(k))/sc
               w1(i,j,k)= we(t-dt,rnu,x(i),y(j),z(k))/sc
               p(i,j,k)= pe(rnu,t/2.0,x(i),y(j),z(k))/sc
               p1(i,j,k)=pe(rnu,(t-dt)/2.0,x(i),y(j),z(k))/sc
               u(i,j,k)= ue(t,rnu,x(i),y(j),z(k))/sc
               v(i,j,k)= ve(t,rnu,x(i),y(j),z(k))/sc
               w(i,j,k)= we(t,rnu,x(i),y(j),z(k))/sc
            end do
         end do
      end do
               ifsetpressure=.false.
               ifzero=.false.
               else
      do k=0,l
         do j=0,n
            do i=0,m
               u(i,j,k)=0.0d0
               v(i,j,k)=0.0d0
               w(i,j,k)=0.0d0
               p(i,j,k)=0.0d0
               p1(i,j,k)=0.0d0
               u1(i,j,k)=0.0d0
               v1(i,j,k)=0.0d0
               w1(i,j,k)=0.0d0
            end do
         end do
      end do
               ifsetpressure=.true.
               ifzero = .true.
               endif  
      !       u1=0.d0;v1=0.d0;w1=0.d0;u=0.d0;v=0.d0;w=0.d0;p=0.d0;p1=0.d0

      !initialize time marching

      !       write(6,*) 'step, time, energy'
      time = dtime(tmp)
      ki = 0
      ! suj is average of uj --CQ

      !       initialize level set
       if(numbubble .eq. 1) then
       call phiset(m,n,l,r0,x0,y0,z0,x,y,z,phi)
       else if(numbubble .eq. 2) then
       call phiset2(m,n,l,x,y,z,phi)
       else
       call phiset3(m,n,l,x,y,z,phi)
       endif

   !  do k=0,l
   !     do j=0,n
   !        do i=0,m
   !        phi0(i,j,k)=-phi(i,j,k)
   !        end do
   !     end do
   !  end do
      ! restart initialization
      ! now u1 and p1 are not used
      ! or can get them from the former step
      ! read the pause points
      restart = 0
      if (restart == 1) then
         ki = 1
         open(unit=1,file='u.restart')
         read(1,*) u
         close(1)
         open(unit=1,file='v.restart')
         read(1,*) v
         close(1)
         open(unit=1,file='w.restart')
         read(1,*) w
         close(1)
         open(unit=1,file='p.restart')
         read(1,*) p
         close(1)
         open(unit=1,file='u1.restart')
         read(1,*) u1
         close(1)
         open(unit=1,file='v1.restart')
         read(1,*) v1
         close(1)
         open(unit=1,file='w1.restart')
         read(1,*) w1
         close(1)
         open(unit=1,file='p1.restart')
         read(1,*) p1
         close(1)
         open(unit=1,file='pi.restart')
         read(1,*) phi
         close(1)
         ifsetpressure=.false.
         ifzero=.false.
      end if
      
      !starting time matching
      !7777   do 1111 ki=1,kmax
      7777 do while (t < tfinal)

        call nsetimer_start(7)
        call nsetimer_start(10)
         ! detect domains
         al = 0.0d0*h
      
           call divide(al,m,n,l,phi,kzone,nzone)


   if(numbubble .eq. 3) then
      do j=1,3
      ctfactor(j)=0.0d0
      enddo     
!   do j=1,3
!   izonex=nint((zonexyz(j,1)-xa)/hx)
!   izoney=nint((zonexyz(j,2)-ya)/hy)
!   izonez=nint((zonexyz(j,3)-za)/hz)
!   tolphi=phi(izonex,izoney,izonez)
!       do k1=-1,1
!           do j1=-1,1
!             do i1=-1,1
!         if(phi(izonex+i1,izoney+j1,izonez+k1) .lt. tolphi) then
!         zonexyz(j,1)=(izonex+i1)*hx+xa  
!         zonexyz(j,2)=(izonex+j1)*hy+ya  
!         zonexyz(j,3)=(izonex+k1)*hz+za  
!         tolphi=phi(izonex+i1,izoney+j1,izonez+k1)
!         endif
!                enddo
!             enddo
!         enddo
!    enddo
          
            
     do j=1,3
     izonex=nint((zonexyz(j,1)-xa)/hx)
     izoney=nint((zonexyz(j,2)-ya)/hy)
     izonez=nint((zonexyz(j,3)-za)/hz)
     k=kzone(izonex,izoney,izonez)/100000
     ctfactor(k)=ctfactor(k)+1.0d0
     if(master) then
     write(6,*) k,ctfactor(k)
     endif
     enddo
     endif
    !  stop       
      !   if (master) print *,'ki',ki,'nzone',nzone

         ! compute volumes
         al = 2.0d0*h
       do ii = 1, nzone
            call setbnd(al,m,n,l,ii,phi,kzone)
            aret(ii) = vol(phi,kzone,ii,al,hx,hy,hz,m,n,l)
            if(master) then
            write(6,*),'volume,step',ii,aret(ii)*125
            endif
        
       end do
         
        
         ! detect projection points
         call indexing2(m,n,l,n2,index2,phi)
         ! information on the projection points
         n3 = 3*n2
         allocate (cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2),cinfo(np*n2)&
                   ,fnt(n3),bf(n3),fvec(n3))
         call curvinf2(m,n,l,n2,np,index2,xa,xb,ya,yb,za,zb,h,hx,hy,hz &
              ,x,y,z,phi,cinfo,kzone)
         ! detect and compute the intersection point of grid and interface
         call bound_info(m,n,l,np,n2,xa,xb,ya,yb,za,zb,hx,hy,hz &
              ,x,y,z,phi,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,kzone)  

     
         ! initialize the vectors for jump condition

         !       call fset(m,n,n2,np,index2,t-dt,rnu,phi,cinfo,ft1,fn1)
         !       call fset(m,n,n2,np,index2,t,rnu,phi,cinfo,ft,fn)
         !       call fset(m,n,n2,np,index2,t+dt,rnu,phi,cinfo,ft2,fn2)

         !       if(ki/3*3.eq.ki) write(6,*)  ki, t, eng


         !       call projout(ki,m,n,l,np,n2,a,b,c,d,h,x,y,phi,
         !    1          indexx,indexy,cinfox,cinfoy,t)
                kll = ki + 9
         !       write(str,'(i6)') kll
         !       str = adjustl(str)
         !       write(fname,'(a)') "uv."//str
         !       open(unit=9,file=fname)
         !       do j=0,n
         !          do i=0,m
         !             if (phi(i,j) .ge. 0.d0) then
         !                write(9,*) u(i,j),v(i,j)
         !             else
         !                write(9,*) 0.d0,0.d0
         !             end if
         !          end do
         !       end do
         !       close(9)
         !       stop
         
         ! compute velocity, pressure
         mm=n3
        call nsetimer_stop(10)
        call nsetimer_stop(7)
 
         call one_step(m,n,l,np,n2,n3,mm,ki,indexx,indexy,indexz &
               ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p &
               ,fw,gw,hw,u1,v1,w1,p1,phi &
               ,cinfox,cinfoy,cinfoz,index2,cinfo &
               ,fnt,bf,fvec,kzone,nzone,aret)
         call mpi_barrier(mpi_comm_world,ierr)
         
         if(nse_debug .le. 8 .and. nse_debug .ne.0) then
         if(t<1.0d0/10.0d0) then
         deallocate(cinfox,cinfoy,cinfoz,cinfo,fnt,bf,fvec)
         dt=t_standard
         t=t+dt
         go to 600
         else
         deallocate(cinfox,cinfoy,cinfoz,cinfo,fnt,bf,fvec)
         return
         endif
         endif
        call nsetimer_start(7)
        call nsetimer_start(10)
         ! !update the results of velocity and pressure
         ! do k=0,l
         !    do j=0,n
         !       do i=0,m
         !          u1(i,j,k)=u(i,j,k)
         !          v1(i,j,k)=v(i,j,k)
         !          w1(i,j,k)=w(i,j,k)  
         !          p1(i,j,k)=p(i,j,k)
         !       end do
         !    end do
         ! end do

         !per-step printing
         if(master) then
        write(6,*)'dt, time, eng,iter ', dt*30.0d0,t*30.0d0,eng,kll
         endif
         ki = ki + 1

         ! print out interface --CQ
         if (master) then
            if( (ki-1) == int((ki-1)/kout)*kout) then
               call projout(ki,m,n,l,np,n2,xa,xb,ya,yb,za,zb,h,x,y,z,phi &
                     ,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,t)
               !          call phiout(m,n,phi)
            end if
         end if

         deallocate(cinfox,cinfoy,cinfoz,cinfo,fnt,bf,fvec)
         ! print the normal direction of the interface
         kll = ki + 9
         str = adjustl(str)

       !  if (master) then
        !    write(fname ,'(a)') "gfi."//str
        !    write(fname1,'(a)') "fix."//str
        !    write(fname2,'(a)') "fiy."//str
        !    write(fname3,'(a)') "fiz."//str
        !    open(unit=9 ,file=fname )
        !    open(unit=10,file=fname1)
        !    open(unit=11,file=fname2)
        !    open(unit=12,file=fname3)
       !  end if
         do k=0,l
            do j=0,n
               do i=0,m
                  call phiinf1(m,n,l,i,j,k,hx,hy,hz,phi,phx,phy,phz)
                  dle = sqrt(phx*phx + phy*phy + phz*phz)
                  coc1 = phx/(dle + 1.0d-10)
                  coc2 = phy/(dle + 1.0d-10)
                  coc3 = phz/(dle + 1.0d-10)
                  un(i,j,k) = u(i,j,k)*coc1 + v(i,j,k)*coc2 &
                        + w(i,j,k)*coc3
                 ! if (master) then
         !            write(9 ,'(f12.6)') dle
         !            write(10,'(f12.6)') coc1
         !            write(11,'(f12.6)') coc2
         !            write(12,'(f12.6)') coc3
                !  end if
               end do
            end do
         end do
       !  if (master) then
         !   close(9 )
         !   close(10)
         !   close(11)
         !   close(12)
        ! end if
!      do k=0,l
!         do j=0,n
!            do i=0,m
!            phi(i,j,k)=phi0(i,j,k)
!            end do
!         end do
!      end do

         !update level set function
         !vhaja:calculate velocity on normal direction
         !advance:recalculate of level set
         !reinit:reinitialization of level set
         width = 20.0*h
         width2 = 25.0*h
         call vhaja(m,n,l,hx,hy,hz,width,phi,un,haja,snorm)
         call advance(m,n,l,dt,phi,phi,haja)
        
         al=0.0d0*h
         call divide(al,m,n,l,phi,kzone,nzone)
         al = 2.0d0*h
         do ii = 1, nzone
            call setbnd(al,m,n,l,ii,phi,kzone)
            aret(ii) = vol(phi,kzone,ii,al,hx,hy,hz,m,n,l)
            if(master) then
            write(6,*),'volume',ii,aret(ii)
            endif
         if(aret(ii)<0.2) then  
          do k=1,l-1
             do j=1,n-1
                 do i=1,m-1
                  if(kzone(i,j,k)==ii*100000) then 
                  kdm=ii 
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                        ,phi,x(i),y(j),z(k),x,y,z,u,rec_value,kdm,kzone)
                  u(i,j,k)=rec_value    
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                        ,phi,x(i),y(j),z(k),x,y,z,v,rec_value,kdm,kzone)
                  v(i,j,k)=rec_value    
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                        ,phi,x(i),y(j),z(k),x,y,z,w,rec_value,kdm,kzone)
                  w(i,j,k)=rec_value    
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                        ,phi,x(i),y(j),z(k),x,y,z,p,rec_value,kdm,kzone)
                  p(i,j,k)=rec_value    
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                        ,phi,x(i),y(j),z(k),x,y,z,u1,rec_value,kdm,kzone)
    !            u1(i,j,k)=rec_value    
    !            call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
    !                  ,phi,x(i),y(j),z(k),x,y,z,v1,rec_value,kdm,kzone)
    !            v1(i,j,k)=rec_value    
    !            call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
    !                  ,phi,x(i),y(j),z(k),x,y,z,w1,rec_value,kdm,kzone)
    !            w1(i,j,k)=rec_value    
    !            call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
    !                  ,phi,x(i),y(j),z(k),x,y,z,p1,rec_value,kdm,kzone)
    !            p1(i,j,k)=rec_value    
                call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i,j,k &
                      ,phi,x(i),y(j),z(k),x,y,z,phi,rec_value,kdm,kzone)
                  phi(i,j,k)=rec_value
                  endif  
                 end do
              end do
           end do
         endif   
         end do

         if( (ki) == int((ki)/kout1)*kout1) then
            al = 0.0
            num = 1
            do ii = 1,5
               call reinit(m,n,l,al,hx,hy,hz,1,1,phi,work,isign,width2,num)
            end do
         end if
 !    do k=0,l
 !       do j=0,n
 !          do i=0,m
 !          phi0(i,j,k)=phi(i,j,k)
 !          end do
 !       end do
 !    end do

         !calculate  energy

         eng=0.0d0
         umax = 0.0d0
         do k=1,l-1
            do j=1,n-1
               do i=1,m-1
                  if( phi(i,j,k) > 0.0) then
                     unorm = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) &
                           + w(i,j,k)*w(i,j,k)
                     eng = eng + unorm
                     unorm = sqrt(unorm)
                     if(umax < unorm) umax = unorm
                  end if
               end do
            end do
         end do
         eng=sqrt(eng*hx*hy*hz)

         !adjusting time step
         t=t+dt
         dt0=dt
         dt = min(t_standard,0.5d0*h/umax)
         ratio_t=dt/dt0

        call nsetimer_stop(10)
        call nsetimer_stop(7)
        600 continue
          

         if(master) then
         call nsetimer_summary(10)
         endif
         ! print out normal velocity --CQ
         ! for plotting --CQ
      !   if (master) then

     !      write(fname,'(a)') "un."//str
     !      open(unit=9,file=fname)
     !      do k=0,l
     !         do j=0,n
     !            do i=0,m
     !               if (phi(i,j,k) > 0.d0) then
     !                  write(9,'(f12.6)') un(i,j,k)
     !               else
     !                  ! write(9,*) 0.d0,0.d0
     !                  write(9,'(f12.6)') un(i,j,k)
     !               end if
     !            end do
     !         end do
     !      end do
     !      close(9)
           if(master) then
           !print velocity x
           write(fname,'(a,I3.2)') "u.",kll
           open(unit=9,file=fname)
           k=int((l+1)/2)
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') u(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') u(i,j,k)
                    end if
                 end do
              end do
     !      end do
           close(9)

           !print velocity y
           write(fname,'(a,I3.2)') "v.",kll
           open(unit=9,file=fname)
           k=int((l+1)/2)
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') v(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') v(i,j,k)
                    end if
                 end do
              end do
           close(9)
          
           if(restart==0) then
           write(fname,'(a)') "u.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') u(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') u(i,j,k)
                    end if
                 end do
              end do
           end do
           close(9)

           !print velocity y
           write(fname,'(a)') "v.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') v(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') v(i,j,k)
                    end if
                 end do
              end do
           end do
           close(9)

           !print velocity z
          write(fname,'(a)') "w.restart"
          open(unit=9,file=fname)
          do k=0,l
             do j=0,n
                do i=0,m
                   if (phi(i,j,k) > 0.d0) then
                      write(9,'(f12.6)') w(i,j,k)
                   else
                      !                write(9,*) 0.d0,0.d0
                      write(9,'(f12.6)') w(i,j,k)
                   end if
                end do
             end do
          end do
          close(9)
 
           write(fname,'(a)') "u1.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') u1(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') u1(i,j,k)
                    end if
                 end do
              end do
           end do
           close(9)

           !print velocity y
           write(fname,'(a)') "v1.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    if (phi(i,j,k) > 0.d0) then
                       write(9,'(f12.6)') v1(i,j,k)
                    else
                       !                write(9,*) 0.d0,0.d0
                       write(9,'(f12.6)') v1(i,j,k)
                    end if
                 end do
              end do
           end do
           close(9)

           !print velocity z
          write(fname,'(a)') "w1.restart"
          open(unit=9,file=fname)
          do k=0,l
             do j=0,n
                do i=0,m
                   if (phi(i,j,k) > 0.d0) then
                      write(9,'(f12.6)') w1(i,j,k)
                   else
                      !                write(9,*) 0.d0,0.d0
                      write(9,'(f12.6)') w1(i,j,k)
                   end if
                end do
             end do
          end do
          close(9)


           ! print out pressure --CQ
           ! for plotting --CQ
           write(fname,'(a)') "p.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    write(9,'(f12.6)') p(i,j,k)
                 end do
              end do
           end do
           close(9)

           ! print out pressure --CQ
           ! for plotting --CQ
           write(fname,'(a)') "p1.restart"
           open(unit=9,file=fname)
           do k=0,l
              do j=0,n
                 do i=0,m
                    write(9,'(f12.6)') p1(i,j,k)
                 end do
              end do
           end do
           close(9)

            
          write(fname,'(a)') "pi.restart"
          open(unit=9,file=fname)
          do k=0,l
             do j=0,n
                do i=0,m
                   write(9,'(f12.6)') phi(i,j,k)
                end do
             end do
          end do
          close(9)
     !     endif

           write(fname,'(a)') "time"
           open(unit=9,file=fname)
           write(9,*)'dt, time, eng,iter ', dt*30.0d0,t*30.0d0,eng,kll
           close(9)
         endif !(restart)
        end if  ! (master)
         ! To see if the normal velocity are the same for all directions for
         ! single sphere
            if (master) then
            nb = 0
            tol_pt = 0.d0
            tol_pt2 = 0.d0
            unmax = -1000.d0
            unmin = 1000.d0
            do k=1,l-1
               do j=1,n-1
                  do i=1,m-1
                     n22 = index2(i,j,k)
                     if(n22 < 0) then
                        nb = nb + 1
                        pt = haja(i,j,k)
                        tol_pt = tol_pt + pt
                        tol_pt2 = tol_pt2 + pt*pt
                        if (pt > unmax) unmax=pt
                        if (pt < unmin) unmin=pt
                        !              unn(nb) = pt
                        !              print *,'unn',nb,pt
                     end if
                  end do
               end do
            end do
            pt = tol_pt / nb
            pt_std = sqrt(tol_pt2/nb-pt*pt)
            !       flux = r0*pt
            !       r0 = r0 + pt*dt
            write(6,*),'ub=',pt,'nbound',nb
            write(6,*),'std',pt_std,'max',unmax,'min',unmin
            !       print *,'r0=',r0,
            !    1   'ub0=',(-2.d0+sqrt(4.d0-2*r0*r0+2*r0))/r0
            !    1   'ub0=',(2.d0-sqrt(4.d0-2*r0*r0+2*r0))/r0
            !    1   'ub0=',(2.0*rnu-sqrt(4.0*rnu*rnu-2.0*(1.0-r0)))/r0
            end if  ! (master)
          !       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !       stop
      end do  ! while (t < tfinal)
     ! 1111 !UNUSED:   continue

#ifdef MPI
 !     call mpi_finalize(ierr)
#endif

      ! check error in u,v,w,p
      ! u,v,w: velocity obtained on three directions
      ! ue,ve,we: analytical velocity obtained on three directions
      ! p:pressure obtained
      ! pe: analytical pressure

   !  time = dtime(tmp)
   !  e1 = 0.0
   !  e2 = 0.0
   !  e3 = 0.0
   !  e4 = 0.0
   !  do k=0,l
   !     do j=0,n
   !        do i=0,m
   !           ui(i,j,k) = 0.0
   !           vi(i,j,k) = 0.0
   !           wi(i,j,k) =0.0      
   !           ! only outside --CQ
   !           if( phi(i,j,k) > 0.0) then
   !              et = abs(u(i,j,k)-ue(t,rnu,x(i),y(j),z(k)))
   !              if(et > e1) e1 = et
   !              et = abs(v(i,j,k)-ve(t,rnu,x(i),y(j),z(k)))
   !              if(et > e2) e2 = et
   !              et = abs(w(i,j,k)-we(t,rnu,x(i),y(j),z(k)))
   !              if(et > e3) then
   !                 e3 = et
   !                 i1 = i
   !                 j1 = j
   !                 k1 = k
   !              end if
   !              et = abs(p(i,j,k)-pe(rnu,t,x(i),y(j),z(k)))
   !              if(et > e4) e4 = et
   !           end if
   !        end do
   !     end do
   !  end do

      !        write(*,*)'error in u,v,w,p ', m,n,l,e1,e2,e3,i1,j1,k1,
      !    1        e1+e2+e3,e4
      !        write(*,*) 'time for iim', time

      100 format(513e16.6)

   end subroutine

   !# ++++printing++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiout here]
   subroutine phiout(m,n,l,phi)
      implicit none
      _REAL_ phi(0:m,0:n,0:l)
      integer m,n,l,i,j,k

      open(5,file='do.m',status='unknown')
      do k=0,l
         do j=0,n
            write(5,100)(phi(i,j,k),i=0,m)
         end do
      end do
      close(5)

      100 format(521e24.16)

      return
   end subroutine

   !+++check++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine checku here]
   subroutine checku(m,n,l,t,x,y,z,phi,u,e,i0,j0,k0,rnu)
      implicit none
      _REAL_ u(0:m,0:n,0:l),x(0:m),y(0:n),z(0:l) ,phi(0:m,0:n,0:l)
      integer m,n,l,i0,j0,k0,i,j,k
      _REAL_ t,e,rnu,et 

      
      et=0
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               ! only outside --CQ
               
               if(phi(i,j,k) > 0.0) then
                  et = -(ue(t,rnu,x(i),y(j),z(k))-u(i,j,k))
                  if( abs(et) > abs(e)) then
                     e = et
                     i0 = i
                     j0 = j
                     k0 = k
                  end if
               end if
            end do
         end do
      end do

      return
   end subroutine

   !++check++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine checkv here]
   subroutine checkv(m,n,l,t,x,y,z,phi,u,e,i0,j0,k0,rnu)
      implicit none
      _REAL_ u(0:m,0:n,0:l),x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      integer m,n,l,i0,j0,k0,i,j,k
      _REAL_ e,et,t,rnu,r0

      e = 0.0
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               ! only outside --CQ
               if(phi(i,j,k) > 0.0) then
                  et = -(ve(t,rnu,x(i),y(j),z(k))-u(i,j,k))
                  if( abs(et) > abs(e)) then
                     e = et
                     i0 = i
                     j0 = j
                     k0 = k
                  end if
               end if
            end do
         end do
      end do

      return
   end subroutine

   !++++++check+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine checkw here]
   subroutine checkw(m,n,l,t,x,y,z,phi,u,e,i0,j0,k0,rnu)
      implicit none
      _REAL_ u(0:m,0:n,0:l),x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      _REAL_ e,et,t,rnu,r0
      integer m,n,l,i0,j0,k0,i,j,k

      e = 0.0
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               ! only outside --CQ
               if(phi(i,j,k) > 0.0) then
                  et = -(we(t,rnu,x(i),y(j),z(k))-u(i,j,k))
                  if( abs(et) > abs(e)) then
                     e = et
                     i0 = i
                     j0 = j
                     k0 = k
                  end if
               end if
            end do
         end do
      end do

      return
   end subroutine

   !# ++++++printing++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine iphiout here]
   subroutine iphiout(m,n,l,iphi)
      implicit none
      integer iphi(0:m,0:n,0:l)
      integer m,n,l,i,j,k

      open(5,file='do1.m',status='unknown')
      do k=0,l
         do j=0,n
            write(5,110)(iphi(i,j,k),i=0,m)
         end do
      end do
      close(5)

      110 format(521i6)

      return
   end subroutine

   !# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !*************************************************************************
   !             SUBROUTINE PRJOUT
   !******printing*******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine projout here]
   subroutine projout(kk,m,n,l,np,n2,xa,xb,ya,yb,za,zb,h,x,y,z &
         ,phi,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,t)

      implicit none

      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
      _REAL_ cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2)
      integer kk,m,n,l,n2,np,i,j,k,kll,ninsc,nx
      _REAL_  xa,xb,ya,yb,za,zb,h,t,f3,fc,fl,rmax,rmin,rr,tot_rr
      _REAL_  x1,x2,x3,xi,yj,zk

 
      kll = kk + 9
      open(kll)
      write(kll,*) 'time',t
      tot_rr = 0.d0
      rmax = 0.d0
      rmin = 100.d0
      ninsc = 0

      do 50 k=2,l-1
         do 50 j=2,n-1
            do 50 i=2,m-1
               !.................. between x_{i-1} and x_i ...........................

               x1 = x(i) - h
               x2 = x(i)
               x3 = x(i) + h
               fl = phi(i-1,j,k)
               fc = phi(i,j,k)
               f3 = phi(i+1,j,k)
               !......                                         !# Irregular grid?
               ! interface --CQ
               if( fl > 0.0 .and. fc <= 0.0) then
                  nx = indexx(i,j,k)
                  xi = cinfox(nx*np-np+1)
                  ninsc = ninsc+1
                  rr = sqrt(xi*xi+y(j)*y(j)+z(k)*z(k))
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)xi,y(j),z(k)
               end if

               if( f3 > 0.0 .and. fc <= 0.0) then
                  if(phi(i-1,j,k) > 0.0) then
                     nx = indexx(i,j,k) + 1
                  else
                     nx = indexx(i,j,k)
                  end if
                  xi = cinfox(nx*np-np+1)
                  ninsc = ninsc+1
                  rr = sqrt(xi*xi+y(j)*y(j)+z(k)*z(k))
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)xi,y(j),z(k)
               end if

      50 continue

      !--- Correct in y-direction ----------------------------------------

      do 60 k=2,l-1
         do 60 j=2,n-1
            do 60 i=2,m-1
               !.................. between y_{j-1} and y_j ...........................

               x1 = y(j) - h
               x2 = y(j)
               x3 = y(j) + h
               fl = phi(i,j-1,k)
               fc = phi(i,j,k)
               f3 = phi(i,j+1,k)
               !......                                         !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  nx = indexy(i,j,k)
                  yj = cinfoy(nx*np-np+2)
                  ninsc = ninsc+1
                  rr = sqrt(x(i)*x(i)+yj*yj+z(k)*z(k))
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)x(i),yj,z(k)
                  !          write(*,'(a,3i6,2f12.6)') 'proj1',i,j,k,fl,fc
               end if

               if( f3 > 0.0 .and. fc <= 0.0) then
                  if(phi(i,j-1,k) > 0.0) then
                     nx = indexy(i,j,k) + 1
                  else
                     nx = indexy(i,j,k)
                  end if
                  yj = cinfoy(nx*np-np+2)
                  ninsc = ninsc+1
                  rr = sqrt(x(i)*x(i)+yj*yj+z(k)*z(k))
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)x(i),yj,z(k)
                  !          write(*,'(a,3i6,2f12.6)') 'proj2',i,j,k,f3,fc
               end if

      60 continue

      !--- Correct in z-direction ----------------------------------------

      do 70 k=2,l-1
         do 70 j=2,n-1
            do 70 i=2,m-1
               !.................. between z_{j-1} and z_j ...........................

               x1 = z(k) - h
               x2 = z(k)
               x3 = z(k) + h
               fl = phi(i,j,k-1)
               fc = phi(i,j,k)
               f3 = phi(i,j,k+1)
               !......                                         !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  nx = indexz(i,j,k)
                  zk = cinfoz(nx*np-np+3)
                  ninsc = ninsc+1
                  rr = sqrt(x(i)*x(i)+y(j)*y(j)+zk*zk)
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)x(i),y(j),zk
               end if

               if( f3 > 0.0 .and. fc <= 0.0) then
                  if(phi(i,j,k-1) > 0.0) then
                     nx = indexz(i,j,k) + 1
                  else
                     nx = indexz(i,j,k)
                  end if
                  zk = cinfoz(nx*np-np+3)
                  ninsc = ninsc+1
                  rr = sqrt(x(i)*x(i)+y(j)*y(j)+zk*zk)
                  if ( rr > rmax ) rmax=rr
                  if ( rr < rmin ) rmin=rr
                  tot_rr = tot_rr+rr
                  write(kll,141)x(i),y(j),zk
               end if

      70 continue

      r0 = tot_rr/ninsc
      write(6,*),'r0=',r0,rmax,rmin

      close(kll)

      141 format(3f12.6)
      return
   end subroutine

   !# +domain detection+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine divide here]
   subroutine divide(al,m,n,l,phi,kzone,nzone)

      implicit none
      _REAL_ phi(0:m,0:n,0:l)
      integer kzone(0:m,0:n,0:l)
      integer i,j,k,nzone,m,n,l,i0,j0,k0
      _REAL_ al
  
      nzone = 0
      kzone = -1
     
      do k = 0, l
         do j = 0, n
            do i = 0, m
               if ( kzone(i,j,k) .ne. -1 ) cycle
               if ( phi(i,j,k) > al ) then
                  kzone(i,j,k) = 0
               else
                  i0=i
                  j0=j
                  k0=k  
                  nzone = nzone + 1
                  call walk(al,i0,j0,k0,m,n,l,phi,kzone,nzone)
               end if
            end do
         end do
      end do

   end subroutine

   !# +dfs depth for search algorithm+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   recursive subroutine walk(al,i,j,k,m,n,l,phi,kzone,nzone)

   implicit none
   _REAL_ phi(0:m,0:n,0:l)
   integer kzone(0:m,0:n,0:l)
   integer i,j,k,m,n,l,nzone
   _REAL_ al   

   kzone(i,j,k) = nzone*100000

   if ( kzone(i+1,j,k) == -1 .and. phi(i+1,j,k) <= al ) &
         call walk(al,i+1,j,k,m,n,l,phi,kzone,nzone)
   if ( kzone(i-1,j,k) == -1 .and. phi(i-1,j,k) <= al ) &
         call walk(al,i-1,j,k,m,n,l,phi,kzone,nzone)
   if ( kzone(i,j+1,k) == -1 .and. phi(i,j+1,k) <= al ) &
         call walk(al,i,j+1,k,m,n,l,phi,kzone,nzone)
   if ( kzone(i,j-1,k) == -1 .and. phi(i,j-1,k) <= al ) &
         call walk(al,i,j-1,k,m,n,l,phi,kzone,nzone)
   if ( kzone(i,j,k+1) == -1 .and. phi(i,j,k+1) <= al ) &
         call walk(al,i,j,k+1,m,n,l,phi,kzone,nzone)
   if ( kzone(i,j,k-1) == -1 .and. phi(i,j,k-1) <= al ) &
         call walk(al,i,j,k-1,m,n,l,phi,kzone,nzone)

end subroutine
!FIXME:                nse.f, line  877: bad indentation level at end of subroutine divide

   !# +set boundary+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setbnd here]
   subroutine setbnd(al,m,n,l,idx,phi,kzone)

      implicit none
      _REAL_ phi(0:m,0:n,0:l)
      integer kzone(0:m,0:n,0:l)
      integer m,n,l,idx,i,j,k,jdx,kk
      _REAL_ al

      jdx = idx*100000
      do k = 0, l
         do j = 0, n
            do i = 0, m
               if ( kzone(i,j,k) /= jdx ) cycle
               kk = kzone(i+1,j,k)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i+1,j,k) = kzone(i+1,j,k) + idx
               kk = kzone(i-1,j,k)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i-1,j,k) = kzone(i-1,j,k) + idx
               kk = kzone(i,j+1,k)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i,j+1,k) = kzone(i,j+1,k) + idx
               kk = kzone(i,j-1,k)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i,j-1,k) = kzone(i,j-1,k) + idx
               kk = kzone(i,j,k+1)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i,j,k+1) = kzone(i,j,k+1) + idx
               kk = kzone(i,j,k-1)
               if ( kk < jdx .and. match(kk,idx) == 0 ) &
                     kzone(i,j,k-1) = kzone(i,j,k-1) + idx
            end do
         end do
      end do

   end subroutine

   !# intersection of domains on the interface++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of integer function match here]
   integer function match(mm,nn)

      integer mm,nn
      integer i,n0,n1,m1,m2,m3

      n0 = 1
      do i = 1, nn-1
         n0 = n0*2
      end do

      n1 = n0*2
      m1 = mod(mm,n1)
      if ( n0 == 1 ) then
         if ( m1 == n0 ) then
            match = 1
            return
         else
            match = 0
            return
         end if
      else
         m2 = mod(m1,n0)
         m3 = m1 - m2
         if ( m3 == n0 ) then
            match = 1
            return
         else
            match = 0
            return
         end if
      end if

   end function


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine advance here]
   subroutine advance(m,n,l,dt,phi,phi1,haja)
      implicit none

      _REAL_ phi(0:m,0:n,0:l),phi1(0:m,0:n,0:l),haja(0:m,0:n,0:l)
      integer i,j,k,i0,j0,k0,m,n,l
      _REAL_ dt

      do 10 k=0,l
         do 10 j=0,n
            do 10 i=0,m
               phi1(i,j,k) = phi(i,j,k) + dt*haja(i,j,k)
      10 continue

      return
   end subroutine

   !*************************************************************************
   !             SUBROUTINE VHAJA
   !*************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine vhaja here]
   subroutine vhaja(nx,ny,nz,dx,dy,dz,width,fi,v,haja,snorm)
      implicit none

      _REAL_ fi(nx+1,ny+1,nz+1),haja(nx+1,ny+1,nz+1) &
            ,v(nx+1,ny+1,nz+1),snorm(nx,ny,nz,2)
      integer nx,ny,nz   
      _REAL_ dx,dy,dz,width,d1x,d1y,d1z,d2x,d2y,dlx,dly,dlz,drx,dry,drz
      _REAL_ dx2,dy2,dz2,ghost,d2z
      integer i,j,k,l,no 

     l = 1
      dx2=dx*dx
      dy2=dy*dy
      dz2=dz*dz

      ! evolution of the level set function
      do k=2,nz
         do j=2,ny
            do i=2,nx

               ! update fi only near the front
               if(dabs(fi(i,j,k)) <= width) then

                  ! upwind scheme by 3rd order WENO

                  ! in x direction

                  if (i == 2) then
                     call weno3(ghost,fi(i-1,j,k),fi(i,j,k),fi(i+1,j,k), &
                           fi(i+2,j,k),dx,i,nx,dlx,drx,d1x,d2x)
                  else if(i == nx) then
                     call weno3(fi(i-2,j,k),fi(i-1,j,k),fi(i,j,k),fi(i+1,j,k), &
                           ghost,dx,i,nx,dlx,drx,d1x,d2x)
                  else
                     call weno3(fi(i-2,j,k),fi(i-1,j,k),fi(i,j,k),fi(i+1,j,k), &
                           fi(i+2,j,k),dx,i,nx,dlx,drx,d1x,d2x)
                  end if

                  ! in y direction

                  if(j == 2) then
                     call weno3(ghost,fi(i,j-1,k),fi(i,j,k),fi(i,j+1,k), &
                           fi(i,j+2,k),dy,j,ny,dly,dry,d1y,d2y)
                  else if(j == ny) then
                     call weno3(fi(i,j-2,k),fi(i,j-1,k),fi(i,j,k),fi(i,j+1,k), &
                           ghost,dy,j,ny,dly,dry,d1y,d2y)
                  else
                     call weno3(fi(i,j-2,k),fi(i,j-1,k),fi(i,j,k),fi(i,j+1,k), &
                           fi(i,j+2,k),dy,j,ny,dly,dry,d1y,d2y)
                  end if

                  ! in z direction

                  if(k == 2) then
                     call weno3(ghost,fi(i,j,k-1),fi(i,j,k),fi(i,j,k+1), &
                           fi(i,j,k+2),dz,k,nz,dlz,drz,d1z,d2z)
                  else if(k == nz) then
                     call weno3(fi(i,j,k-2),fi(i,j,k-1),fi(i,j,k),fi(i,j,k+1), &
                           ghost,dz,k,nz,dlz,drz,d1z,d2z)
                  else
                     call weno3(fi(i,j,k-2),fi(i,j,k-1),fi(i,j,k),fi(i,j,k+1), &
                           fi(i,j,k+2),dz,k,nz,dlz,drz,d1z,d2z)
                  end if

                  ! update fi by upwind scheme for Hamilton Jacobi
                  if(v(i,j,k) > 0.0) then
                     haja(i,j,k)= -v(i,j,k)*dsqrt( &
                           (max(dlx,0.0d0))**2+(min(drx,0.0d0))**2+ &
                           (max(dly,0.0d0))**2+(min(dry,0.0d0))**2+ &
                           (max(dlz,0.0d0))**2+(min(drz,0.0d0))**2)
                  else
                     haja(i,j,k)=-v(i,j,k)*dsqrt( &
                           (min(dlx,0.0d0))**2+(max(drx,0.0d0))**2+ &
                           (min(dly,0.0d0))**2+(max(dry,0.0d0))**2+ &
                           (min(dlz,0.0d0))**2+(max(drz,0.0d0))**2)
                  end if

               end if  ! (dabs(fi(i,j,k)) <= width)
            end do  !  i=2,nx
         end do  !  j=2,ny
      end do  !  k=2,nz

      no = 1

      ! set boundary value by boundary conditions

      ! set boundary by using harmonic Neumann condition
      !       call Neu_BO(fi, no, l, nx, ny)

      ! set boundary by extropolation
      call ex_bo(fi, no, l, nx, ny, nz)

      return
   end subroutine

   !-------------------------------------------------------------------

   ! set boundary value by harmonic Neumann condition
   ! 2D, not called

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine neu_bo here]
   subroutine neu_bo(u, no, l, nx, ny)
      integer no, l, nx, ny
      real*8  u(nx+1,ny+1,no)

      integer i, j
      ! boundaryn
      do i=1,nx+1
         u(i,1,l)=u(i,2,l)
         u(i,ny+1,l)=u(i,ny,l)
      end do

      do j=1,ny+1
         u(1,j,l)=u(2,j,l)
         u(nx+1,j,l)=u(nx,j,l)
      end do

      ! corner points
      u(1,1,l)=u(2,2,l)
      u(1,ny+1,l)=u(2,ny,l)
      u(nx+1,ny+1,l)=u(nx,ny,l)
      u(nx+1,1,l)=u(nx,2,l)

      return
   end subroutine

   !-----------------------------------------------------------------

   ! set boundary by upwind extropolation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ex_bo here]
   subroutine ex_bo(u, no, l, nx, ny, nz)
      integer no, l, nx, ny, nz
      real*8  u(nx+1,ny+1,nz+1,no)
      !       real*8  sgn

      integer i, j, k
      ! side
      do k=1,nz+1
         do j=1,ny+1
            u(1,j,k,l)=u(2,j,k,l)+sgn(u(2,j,k,l)) &
                  *dabs(u(2,j,k,l)-u(3,j,k,l))
            u(nx+1,j,k,l)=u(nx,j,k,l)+sgn(u(nx,j,k,l)) &
                  *dabs(u(nx,j,k,l)-u(nx-1,j,k,l))
         end do
      end do

      do k=1,nz+1
         do i=1,nx+1
            u(i,1,k,l)=u(i,2,k,l)+sgn(u(i,2,k,l)) &
                  *dabs(u(i,2,k,l)-u(i,3,k,l))
            u(i,ny+1,k,l)=u(i,ny,k,l)+sgn(u(i,ny,k,l)) &
                  *dabs(u(i,ny,k,l)-u(i,ny-1,k,l))
         end do
      end do

      do j=1,ny+1
         do i=1,nx+1
            u(i,j,1,l)=u(i,j,2,l)+sgn(u(i,j,2,l)) &
                  *dabs(u(i,j,2,l)-u(i,j,3,l))
            u(i,j,nz+1,l)=u(i,j,nz,l)+sgn(u(i,j,nz,l)) &
                  *dabs(u(i,j,nz,l)-u(i,j,nz-1,l))
         end do
      end do

      ! boundary
      do i=1,nx+1
         u(i,1,1,l)=u(i,2,2,l)+sgn(u(i,2,2,l)) &
               *dabs(u(i,2,2,l)-u(i,3,3,l))
         u(i,ny+1,1,l)=u(i,ny,2,l)+sgn(u(i,ny,2,l)) &
               *dabs(u(i,ny,2,l)-u(i,ny-1,3,l))
         u(i,1,nz+1,l)=u(i,2,nz,l)+sgn(u(i,2,nz,l)) &
               *dabs(u(i,2,nz,l)-u(i,3,nz-1,l))
         u(i,ny+1,nz+1,l)=u(i,ny,nz,l)+sgn(u(i,ny,nz,l)) &
               *dabs(u(i,ny,nz,l)-u(i,ny-1,nz-1,l))
      end do

      do j=1,ny+1
         u(1,j,1,l)=u(2,j,2,l)+sgn(u(2,j,2,l)) &
               *dabs(u(2,j,2,l)-u(3,j,3,l))
         u(nx+1,j,1,l)=u(nx,j,2,l)+sgn(u(nx,j,2,l)) &
               *dabs(u(nx,j,2,l)-u(nx-1,j,3,l))
         u(1,j,nz+1,l)=u(2,j,nz,l)+sgn(u(2,j,nz,l)) &
               *dabs(u(2,j,nz,l)-u(3,j,nz-1,l))
         u(nx+1,j,nz+1,l)=u(nx,j,nz,l)+sgn(u(nx,j,nz,l)) &
               *dabs(u(nx,j,nz,l)-u(nx-1,j,nz-1,l))
      end do

      do k=1,nz+1
         u(1,1,k,l)=u(2,2,k,l)+sgn(u(2,2,k,l)) &
               *dabs(u(2,2,k,l)-u(3,3,k,l))
         u(nx+1,1,k,l)=u(nx,2,k,l)+sgn(u(nx,2,k,l)) &
               *dabs(u(nx,2,k,l)-u(nx-1,3,k,l))
         u(1,ny+1,k,l)=u(2,ny,k,l)+sgn(u(2,ny,k,l)) &
               *dabs(u(2,ny,k,l)-u(3,ny-1,k,l))
         u(nx+1,ny+1,k,l)=u(nx,ny,k,l)+sgn(u(nx,ny,k,l)) &
               *dabs(u(nx,ny,k,l)-u(nx-1,ny-1,k,l))
      end do

      ! corners
      u(1,1,1,l)=u(2,2,2,l)+sgn(u(2,2,2,l)) &
            *dabs(u(2,2,2,l)-u(3,3,3,l))
      u(nx+1,1,1,l)=u(nx,2,2,l)+sgn(u(nx,2,2,l)) &
            *dabs(u(nx,2,2,l)-u(nx-1,3,3,l))
      u(1,ny+1,1,l)=u(2,ny,2,l)+sgn(u(2,ny,2,l)) &
            *dabs(u(2,ny,2,l)-u(3,ny-1,3,l))
      u(nx+1,ny+1,1,l)=u(nx,ny,2,l)+sgn(u(nx,ny,2,l)) &
            *dabs(u(nx,ny,2,l)-u(nx-1,ny-1,3,l))
      u(1,1,nz+1,l)=u(2,2,nz,l)+sgn(u(2,2,nz,l)) &
            *dabs(u(2,2,nz,l)-u(3,3,nz-1,l))
      u(nx+1,1,nz+1,l)=u(nx,2,nz,l)+sgn(u(nx,2,nz,l)) &
            *dabs(u(nx,2,nz,l)-u(nx-1,3,nz-1,l))
      u(1,ny+1,nz+1,l)=u(2,ny,nz,l)+sgn(u(2,ny,nz,l)) &
            *dabs(u(2,ny,nz,l)-u(3,ny-1,nz-1,l))
      u(nx+1,ny+1,nz+1,l)=u(nx,ny,nz,l)+sgn(u(nx,ny,nz,l)) &
            *dabs(u(nx,ny,nz,l)-u(nx-1,ny-1,nz-1,l))

      return
   end subroutine


   !---------------------------------------------------------------

   real*8 function absmin(x, y)
   real*8 x, y

   if (x*y < 0.0d0) then
      absmin = 0.0d0
   else
      if (dabs(x) <= dabs(y)) then
         absmin = x
      else
         absmin = y
      end if
   end if

   return
end function
!FIXME:                nse.f, line 1200: bad indentation level at end of subroutine ex_bo

   !*************************************************************************
   !             SUBROUTINE VHAJA
   !*************************************************************************
   ! 2D, not called




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine arc here]
   subroutine arc(x1,y1,z1,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,x2,y2,z2,cocx,cocy,cocz,sy,sz)

      !# ---------------------------------------------------------------------
      !# - Subroutine arc calculate the signed arc length between two points
      !# - (x1,y1) and (x2,y2) starting from (x1,y1) using Hermite interpolation
      !# -- and Simpson's rule. So it is third order accurate.

      implicit none
      _REAL_ x1,y1,z1,coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ x2,y2,z2,cocx,cocy,cocz,sy,sz,c1,c2,coc1,coc2,coc3,cocb
      _REAL_ cocc,comp,d1,d2,f2y,f2z,fmidy,fmidz,rmag,x3,y3,ymid,ymid1
      _REAL_ z3,zmid,zmid1,coca


      !# ----- Coordinates transform

      !       call gettan(coc11,coc12,coc13,coca1,coca2,coca3
      !    1    ,cocb1,cocb2,cocb3)

      x3 = (x2-x1)*coca1 + (y2-y1)*coca2 + (z2-z1)*coca3
      y3 = (x2-x1)*cocb1 + (y2-y1)*cocb2 + (z2-z1)*cocb3
      z3 = (x2-x1)*cocc1 + (y2-y1)*cocc2 + (z2-z1)*cocc3
      coca = cocx*coca1 + cocy*coca2 + cocz*coca3
      cocb = cocx*cocb1 + cocy*cocb2 + cocz*cocb3
      cocc = cocx*cocc1 + cocy*cocc2 + cocz*cocc3

      comp = coca*x3+cocb*y3+cocc*z3
      coc1 = x3 - comp*coca
      coc2 = y3 - comp*cocb
      coc3 = z3 - comp*cocc
      rmag = 1.d0/dsqrt(coc1*coc1+coc2*coc2+coc3*coc3)
      coc1 = coc1*rmag
      coc2 = coc2*rmag
      coc3 = coc3*rmag
      !       write(54,*) 'test',coc1,coc2,coc3
      !       write(54,*) 'test2',coca*coc1+cocb*coc2+cocc*coc3

      ! for debug --CQ
      !       rmag = 1.d0/dsqrt(x2*x2+y2*y2+z2*z2)
      !       x2 = x2*rmag
      !       y2 = y2*rmag
      !       z2 = z2*rmag
      !       comp = (x2-x1)*x2+(y2-y1)*y2+(z2-z1)*z2
      !       v1 = (x2-x1)-comp*x2
      !       v2 = (y2-y1)-comp*y2
      !       v3 = (z2-z1)-comp*z2
      !       write(54,*) 'test',v1*x2+v2*y2+v3*z2
      !       rmag = 1.d0/dsqrt(v1*v1+v2*v2+v3*v3)
      !       v1 = v1*rmag
      !       v2 = v2*rmag
      !       v3 = v3*rmag
      !       w1 = v1*coca1 + v2*coca2 + v3*coca3
      !       w2 = v1*cocb1 + v2*cocb2 + v3*cocb3
      !       w3 = v1*cocc1 + v2*cocc2 + v3*cocc3
      !       a1=y1*z2-y2*z1
      !       a2=x2*z1-x1*z2
      !       a3=x1*y2-x2*y1
      !       b1 = a1*coca1 + a2*coca2 + a3*coca3
      !       b2 = a1*cocb1 + a2*cocb2 + a3*cocb3
      !       b3 = a1*cocc1 + a2*cocc2 + a3*cocc3
      !       write(54,*) 'test',v1*cocx+v2*cocy+v3*cocz
      !       write(54,*) 'test',v1*x2+v2*y2+v3*z2
      !       write(54,*) 'test',b1*coc1+b2*coc2+b3*coc3
      !       write(54,*) 'test',b1*x3+b2*y3+b3*z3
      !       write(54,'(6f12.6)') w1,w2,w3,coc1,coc2,coc3
      !       stop

      !# ----- Determine the Hermite cubic

      !       if( abs(y3) .le. 1.0d-8 .or. abs(z3) .le. 1.0d-8 ) then
      !          s = dsqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
      !    1       + (z1-z2)*(z1-z2) )
      !       else
      !# the curve is   $ \xi = c \eta^2 + d \eta^3$
      c1 = 3.0d0*x3/(y3*y3) - coc1/(y3*coc2)
      d1 = coc1/(coc2*y3*y3) - 2.0d0*x3/(y3*y3*y3)
      c2 = 3.0d0*x3/(z3*z3) - coc1/(z3*coc3)
      d2 = coc1/(coc3*z3*z3) - 2.0d0*x3/(z3*z3*z3)

      !# --- Using the simpson rule --------------------

      ymid = y3*0.5d0
      ymid1 = ymid*ymid
      zmid = z3*0.5d0
      zmid1 = zmid*zmid
      !         aa = -2.d0*x3/(y3**3)+coc1/(y3*y3*coc2)
      !         bb = 3.d0*x3/(y3*y3)-coc1/(y3*coc2)
      !         cc = -0.5d0*x3
      !         call solvcubic(aa,bb,cc,ymid)
      !         aa = -2.d0*x3/(z3**3)+coc1/(z3*z3*coc3)
      !         bb = 3.d0*x3/(z3*z3)-coc1/(z3*coc3)
      !         cc = -0.5d0*x3
      !         call solvcubic(aa,bb,cc,zmid)
      !         write(54,*) 'mid',ymid,y3,zmid,z3

      !         fmid = dsqrt( 1.d0 +
      !    1      1.d0/(ymid1*(2.d0*c1+3.d0*d1*ymid)*(2.d0*c1+3.d0*d1*ymid)) +
      !    1      1.d0/(zmid1*(2.d0*c2+3.d0*d2*zmid)*(2.d0*c2+3.d0*d2*zmid)) )
      !         f2 = dsqrt( 1.d0 +
      !    1      1.d0/(y3*y3*(2.d0*c1+3.d0*d1*y3)*(2.d0*c1+3.d0*d1*y3)) +
      !    1      1.d0/(z3*z3*(2.d0*c2+3.d0*d2*z3)*(2.d0*c2+3.d0*d2*z3)) )
      !         s= (1.0d0 + 4.0d0*fmid + f2 )*x3/6.0d0
      !         s= (1.d0 + f2)*x3/2.d0

      fmidy = dsqrt( 1.d0 + &
            ymid1*(2.d0*c1+3.d0*d1*ymid)*(2.d0*c1+3.d0*d1*ymid) )
      fmidz = dsqrt( 1.d0 + &
            zmid1*(2.d0*c2+3.d0*d2*zmid)*(2.d0*c2+3.d0*d2*zmid) )
      f2y = dsqrt( 1.d0 + &
            y3*y3*(2.d0*c1+3.d0*d1*y3)*(2.d0*c1+3.d0*d1*y3) )
      f2z = dsqrt( 1.d0 + &
            z3*z3*(2.d0*c2+3.d0*d2*z3)*(2.d0*c2+3.d0*d2*z3) )

      sy = (1.0d0 + 4.0d0*fmidy + f2y )*y3/6.0d0
      sz = (1.0d0 + 4.0d0*fmidz + f2z )*z3/6.0d0
      !       endif

      if( abs(y3) <= 1.0d-8 ) then
         sy = dsqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )
      end if
      if( abs(z3) <= 1.0d-8 ) then
         sz = dsqrt( (x1-x2)*(x1-x2) + (z1-z2)*(z1-z2) )
      end if

      return
   end subroutine



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine solvcubic here]
   subroutine solvcubic(a,b,c,x)

      !       solve a*x^3+b*x^2+c=0

      implicit none
      _REAL_ p,q,a,b,c,s,t,x
      
      

      p = dsqrt(3.d0*(4.d0*a**2*b**3*c+27.d0*a**4*c**2))
      q = (-b**3-13.5d0*a**2*c+1.5d0*p)
      s = abs(q)**(1.d0/3.d0)
      t = sign(s,q)
      x = (-b + b**2/q + q)/a/3.d0
      !       write(54,*) 'cubic',a,b,c,p,q,s,t,x

   end subroutine solvcubic


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of ubroutine bcg here]
  subroutine bcg(m,n,l,np,n2,n3,istart,indexx,indexy,indexz &
         ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p &
         ,fw,gw,hw,u1,v1,w1,p1,phi,cinfox,cinfoy &
         ,cinfoz,index2,cinfo,ft,fn,fm,fj,gj,hj,kzone)

      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,fw(0:m,0:n,0:l),gw(0:m,0:n,0:l),hw(0:m,0:n,0:l) &
            ,ui(0:m,0:n,0:l),vi(0:m,0:n,0:l),wi(0:m,0:n,0:l) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) &
            ,p1(0:m,0:n,0:l),p2(0:m,0:n,0:l),p(0:m,0:n,0:l) &
            ,phi(0:m,0:n,0:l)&
            ,fj(0:m,0:n,0:l),gj(0:m,0:n,0:l),hj(0:m,0:n,0:l)
      integer index2(0:m,0:n,0:l),kzone(0:m,0:n,0:l)
      _REAL_  cinfo(np*n2) &
              ,cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2) 
       integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
       _REAL_ ft(n2),fn(n2),fm(n2)
      integer i,j,k,n3,m,n,l,np,n2,istart
      _REAL_ h,hx1,hx2,hy1,hy2,hz1,hz2,hx,hy,hz,cc,rnu,tt,xa,xb
      _REAL_ ya,yb,za,zb,dt,t,a1,b1,coca1,coca2,coca3,cocb1,cocb2,cocb3
      _REAL_ cocc1,cocc2,cocc3,curv,e,e1,e2,elmbda,f3,fc,fl,fuju,fvjv,fwjw
      integer i0,j0,k0,i1,i2,ierror,ilr,info,j1,j2,k1,k2,ibdcnd,mbdcnd,mdimf
      integer nbdcnd,ndimf,ne,npp,nx,ny,nz,lbdcnd,solvopt,optx,opty,optz,kdm
      _REAL_  pertrb,pnin,pnout,pxp,pyp,pzp,r1,wt,x1,x2,t1,x3,xi,xp,xt,y1,yj,yp
      _REAL_  unout,unin,vnout,vnin,wnout,wnin,uxu,vxv,wxw,uyu,vyv,wyw,uzu,vzv,wzw
      _REAL_  yt,z1,zk,zp,zt,ww(30+10*(m+n+l))
      _REAL_  epsx,epsy,epsz,iv(1:(m-1)*(n-1)*(l-1)),xso((m-1)*(n-1)*(l-1)+2*(m-1)*(n-1)),dummy,accept 
      ! common /radius/r0
      !-------------------------------initialize---------------------
      h = min(hx,hy,hz)
      hx1 = hx*hx
      hx2 = 2.0d0*hx
      hy1 = hy*hy
      hy2 = 2.0d0*hy
      hz1 = hz*hz
      hz2 = 2.0d0*hz
    !  cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      tt = 0.d0
      e=0
      e1=0
      e2=0
      solvopt=1

      call nsetimer_start(1)
      do k=0,l
         do j=0,n
            do i=0,m
               ui(i,j,k) = 0.0d0
               vi(i,j,k) = 0.0d0
               wi(i,j,k) = 0.0d0
               p2(i,j,k)=2.0*p(i,j,k)-p1(i,j,k)
               !                   change algorithm of p
            end do
         end do
      end do   
    call nsetimer_stop(1)
    call nsetimer_start(2)
      if(ifupdate .eqv. .true.) then
      do  k= 1,l-1
         do  j= 1,n-1
            do  i= 1,m-1

               ! interface --CQ
               ! only outside --CQ
               ! index2 only inside --CQ

               if ( phi(i,j,k) > 0.d0) then           
 !# -------------------- Begin x-direction for u*delta(u) ----------------

                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0) then
                     if(phi(i,j,k)*phi(i+1,j,k) > 0.0) then
                     ui(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(u(i+1,j,k)-u(i,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(u1(i+1,j,k)-u1(i,j,k)) /hx
                     vi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(v(i+1,j,k)-v(i,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(v1(i+1,j,k)-v1(i,j,k)) /hx
                     wi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(w(i+1,j,k)-w(i,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(w1(i+1,j,k)-w1(i,j,k)) /hx
                     pxp= (p2(i+1,j,k)-p2(i,j,k))/hx
                     else
                     pxp=0.0d0
                     endif
                     

        !            nx = indexx(i-1,j,k)
        !            if (phi(i-2,j,k) > 0.0) nx = nx+1
        !            x1 = cinfox(nx*np-np+1)
        !            y1 = cinfox(nx*np-np+2)
        !            z1 = cinfox(nx*np-np+3)
        !            coca1 = cinfox(nx*np-np+4)
        !            coca2 = cinfox(nx*np-np+5)
        !            coca3 = cinfox(nx*np-np+6)
        !            cocb1 = cinfox(nx*np-np+7)
        !            cocb2 = cinfox(nx*np-np+8)
        !            cocb3 = cinfox(nx*np-np+9)
        !            cocc1 = cinfox(nx*np-np+10)
        !            cocc2 = cinfox(nx*np-np+11)
        !            cocc3 = cinfox(nx*np-np+12)
        !            kdm=nint(cinfox(nx*np-np+22))
        !            r1 = x1*x1 + y1*y1 + z1*z1
        !            xt = x1/(r1*r1*r1)
        !            yt = y1/(r1*r1*r1)
        !            zt = z1/(r1*r1*r1)
        !            
        !            ! calculate the jump condition of u
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
        !            wt = (unout - unin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            uxu = (u(i,j,k)-u(i-1,j,k)-tt)/hx+a1*wt
        !            ui(i,j,k)=-2.0*u(i,j,k)*uxu


        !            ! calculate the jump condition of u1
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
        !            wt = (unout - unin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            uxu = (u1(i,j,k)-u1(i-1,j,k)-tt)/hx+a1*wt
        !            ui(i,j,k)=ui(i,j,k)+u1(i,j,k)*uxu
        !           
        !            ! calculate the jump condition of v
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
        !            wt = (vnout - vnin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            vxv = (v(i,j,k)-v(i-1,j,k)-tt)/hx+a1*wt
        !            vi(i,j,k)=-2.0*u(i,j,k)*vxv


        !            ! calculate the jump condition of v1
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
        !            wt = (vnout - vnin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            vxv = (v1(i,j,k)-v1(i-1,j,k)-tt)/hx+a1*wt
        !            vi(i,j,k)=vi(i,j,k)+u1(i,j,k)*vxv
        !            
        !           ! calculate the jump condition of w
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
        !            wt = (wnout - wnin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            wxw = (w(i,j,k)-w(i-1,j,k)-tt)/hx+a1*wt
        !            wi(i,j,k)=-2.0*u(i,j,k)*wxw


        !            ! calculate the jump condition of w1
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
        !            wt = (wnout - wnin)*abs(coca1)
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            wxw = (w1(i,j,k)-w1(i-1,j,k)-tt)/hx+a1*wt
        !            wi(i,j,k)=wi(i,j,k)+u1(i,j,k)*wxw

        !            ! calculate the jump condition of pn
        !            call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnout)
        !            call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
        !                  ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
        !                  ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnin,kdm,kzone)
        !            wt = (pnout - pnin)*abs(coca1)
        !            !         print *,'interp4',pnout-pnin
        !            !    1      ,3.d0*cc*cc*(xt*coca1+yt*coca2+zt*coca3)
        !            !         tt = -cc*cc/2.0d0/r1-1.0d0/r1
        !            a1 = (x1-x(i-1))/hx
        !            b1 = 1.d0 - a1
        !            pxp = (p2(i,j,k)-p2(i-1,j,k)-tt)/hx+a1*wt
        !            endif
        !            !         print *, 'pxp-1',i,j,a1,wt,pxp
                     !         print *, 'pxp-1',x1,y1,hx,cc,r1

                  else if(phi(i,j,k)*phi(i+1,j,k) <= 0.0 ) then
                     ui(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(u(i,j,k)-u(i-1,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(u1(i,j,k)-u1(i-1,j,k)) /hx
                     vi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(v(i,j,k)-v(i-1,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(v1(i,j,k)-v1(i-1,j,k)) /hx
                     wi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(w(i,j,k)-w(i-1,j,k))/hx &
                           +ratio_t*u1(i,j,k)*(w1(i,j,k)-w1(i-1,j,k)) /hx

                     pxp= (p2(i,j,k)-p2(i-1,j,k))/hx
                     !         print *, 'pxp+1',i,j,a1,wt,pxp

                  else

                     ui(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(u(i+1,j,k)-u(i-1,j,k))/hx2 &
                           +ratio_t*u1(i,j,k)*(u1(i+1,j,k)-u1(i-1,j,k)) /hx2
                     vi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(v(i+1,j,k)-v(i-1,j,k))/hx2 &
                           +ratio_t*u1(i,j,k)*(v1(i+1,j,k)-v1(i-1,j,k)) /hx2
                     wi(i,j,k)=-(1.0d0+ratio_t)*u(i,j,k)*(w(i+1,j,k)-w(i-1,j,k))/hx2 &
                           +ratio_t*u1(i,j,k)*(w1(i+1,j,k)-w1(i-1,j,k)) /hx2
                     pxp = (p2(i+1,j,k)-p2(i-1,j,k))/hx2

                  end if  ! (phi(i,j,k)*phi(i-1,j,k) <= 0.0)


                  !#  done with x-direction, now in y-direction  ------------------

                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0) then
                    if(phi(i,j,k)*phi(i,j+1,k)>0.0) then

                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(u(i,j+1,k)-u(i,j,k))/hy &
                           +ratio_t*v1(i,j,k)*(u1(i,j+1,k)-u1(i,j,k))/hy
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(v(i,j+1,k)-v(i,j,k))/hy &
                           +ratio_t*v1(i,j,k)*(v1(i,j+1,k)-v1(i,j,k))/hy
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(w(i,j+1,k)-w(i,j,k))/hy &
                           +ratio_t*v1(i,j,k)*(w1(i,j+1,k)-w1(i,j,k))/hy
                     pyp= (p2(i,j+1,k)-p2(i,j,k))/hy
                     else
                     pyp=0.0d0
                     endif
                     

           !         ny = indexy(i,j-1,k)
           !         if (phi(i,j-2,k) > 0.0) ny = ny+1
           !         x1 = cinfoy(ny*np-np+1)
           !         y1 = cinfoy(ny*np-np+2)
           !         z1 = cinfoy(ny*np-np+3)
           !         coca1 = cinfoy(ny*np-np+4)
           !         coca2 = cinfoy(ny*np-np+5)
           !         coca3 = cinfoy(ny*np-np+6)
           !         cocb1 = cinfoy(ny*np-np+7)
           !         cocb2 = cinfoy(ny*np-np+8)
           !         cocb3 = cinfoy(ny*np-np+9)
           !         cocc1 = cinfoy(ny*np-np+10)
           !         cocc2 = cinfoy(ny*np-np+11)
           !         cocc3 = cinfoy(ny*np-np+12)
           !         kdm=nint(cinfoy(ny*np-np+22))
           !         r1 = x1*x1 + y1*y1 + z1*z1
           !         xt = x1/(r1*r1*r1)
           !         yt = y1/(r1*r1*r1)
           !         zt = z1/(r1*r1*r1)
           !         ! calculate the jump condition of u
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
           !         wt = (unout - unin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         uyu = (u(i,j,k)-u(i,j-1,k)-tt)/hy+a1*wt
           !         ui(i,j,k)=ui(i,j,k)-2.0*v(i,j,k)*uyu


           !         ! calculate the jump condition of u1
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
           !         wt = (unout - unin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         uyu = (u1(i,j,k)-u1(i,j-1,k)-tt)/hy+a1*wt
           !         ui(i,j,k)=ui(i,j,k)+v1(i,j,k)*uyu
           !        
           !         ! calculate the jump condition of v
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
           !         wt = (vnout - vnin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         vyv = (v(i,j,k)-v(i,j-1,k)-tt)/hy+a1*wt
           !         vi(i,j,k)=vi(i,j,k)-2.0*v(i,j,k)*vyv


           !         ! calculate the jump condition of v1
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
           !         wt = (vnout - vnin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         vyv = (v1(i,j,k)-v1(i,j-1,k)-tt)/hy+a1*wt
           !         vi(i,j,k)=vi(i,j,k)+v1(i,j,k)*vyv
           !         
           !        ! calculate the jump condition of w
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
           !         wt = (wnout - wnin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         wyw = (w(i,j,k)-w(i,j-1,k)-tt)/hy+a1*wt
           !         wi(i,j,k)=wi(i,j,k)-2.0*v(i,j,k)*wyw


           !         ! calculate the jump condition of w1
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
           !         wt = (wnout - wnin)*abs(coca2)
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         wyw = (w1(i,j,k)-w1(i,j-1,k)-tt)/hy+a1*wt
           !         wi(i,j,k)=wi(i,j,k)+v1(i,j,k)*wyw
           !         !         wt = 3.0d0*cc*cc*(xt*coca1+yt*coca2+zt*coca3)*abs(coca2)
           !         call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnout)
           !         call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
           !               ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
           !               ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnin,kdm,kzone)
           !         wt = (pnout - pnin)*abs(coca2)
           !         !         tt = -cc*cc/2.0d0/r1-1.0d0/r1
           !         a1 = (y1-y(j-1))/hy
           !         b1 = 1.d0 - a1
           !         pyp = (p2(i,j,k)-p2(i,j-1,k)-tt)/hy+a1*wt
           !         endif
                     !         print *, 'pyp-1',i,j,a1,wt,pyp

                  else if(phi(i,j,k)*phi(i,j+1,k) <= 0.0 ) then

                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(u(i,j,k)-u(i,j-1,k))/hy &
                           +ratio_t*v1(i,j,k)*(u1(i,j,k)-u1(i,j-1,k)) /hy
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(v(i,j,k)-v(i,j-1,k))/hy &
                           +ratio_t*v1(i,j,k)*(v1(i,j,k)-v1(i,j-1,k)) /hy
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(w(i,j,k)-w(i,j-1,k))/hy &
                           +ratio_t*v1(i,j,k)*(w1(i,j,k)-w1(i,j-1,k)) /hy

                     pyp= (p2(i,j,k)-p2(i,j-1,k))/hy
                     !         print *, 'pyp+1',i,j,a1,ww,pyp

                     !       else if(phi(i,j)*phi(i-1,j) .le. 0.0 .or.
                     !    1          phi(i,j)*phi(i+1,j) .le. 0.0 ) then
                  else

                     !           ui(i,j,k)= ui(i,j,k)-v(i,j,k)*(u(i,j+1,k)-u(i,j-1,k))/hy2
                     !           vi(i,j,k)= vi(i,j,k)-v(i,j,k)*(v(i,j+1,k)-v(i,j-1,k))/hy2
                     !           wi(i,j,k)= wi(i,j,k)-v(i,j,k)*(w(i,j+1,k)-w(i,j-1,k))/hy2
                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(u(i,j+1,k)-u(i,j-1,k))/hy2 &
                           +ratio_t*v1(i,j,k)*(u1(i,j+1,k)-u1(i,j-1,k)) /hy2
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(v(i,j+1,k)-v(i,j-1,k))/hy2 &
                           +ratio_t*v1(i,j,k)*(v1(i,j+1,k)-v1(i,j-1,k)) /hy2
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*v(i,j,k)*(w(i,j+1,k)-w(i,j-1,k))/hy2 &
                           +ratio_t*v1(i,j,k)*(w1(i,j+1,k)-w1(i,j-1,k)) /hy2
                     pyp = (p2(i,j+1,k)-p2(i,j-1,k))/hy2

                  end if  ! (phi(i,j,k)*phi(i,j-1,k) <= 0.0)
                  !      if (i==16.and.j==14.and.k==9) print *,'y',ui(i,j,k)

                  !#  done with y-direction, now in z-direction  ------------------

                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0) then

                     !         if(phi(i,j,k+1) .gt. 0.0) then
                     !           ui(i,j,k)= ui(i,j,k)-w(i,j,k)*(u(i,j,k+1)-u(i,j,k))/hz
                     !           vi(i,j,k)= vi(i,j,k)-w(i,j,k)*(v(i,j,k+1)-v(i,j,k))/hz
                     !           wi(i,j,k)= wi(i,j,k)-w(i,j,k)*(w(i,j,k+1)-w(i,j,k))/hz
                     if(phi(i,j,k)*phi(i,j,k+1) > 0.0)then
                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(u(i,j,k+1)-u(i,j,k))/hz &
                           +ratio_t*w1(i,j,k)*(u1(i,j,k+1)-u1(i,j,k)) /hz
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(v(i,j,k+1)-v(i,j,k))/hz &
                           +ratio_t*w1(i,j,k)*(v1(i,j,k+1)-v1(i,j,k)) /hz
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(w(i,j,k+1)-w(i,j,k))/hz &
                           +ratio_t*w1(i,j,k)*(w1(i,j,k+1)-w1(i,j,k)) /hz
                     pzp= (p2(i,j,k+1)-p2(i,j,k))/hz
                     else
                     pzp=0.0d0
                     endif
                     !         end if
              !      nz = indexz(i,j,k-1)
              !      if (phi(i,j,k-2) > 0.0) nz = nz+1
              !      x1 = cinfoz(nz*np-np+1)
              !      y1 = cinfoz(nz*np-np+2)
              !      z1 = cinfoz(nz*np-np+3)
              !      coca1 = cinfoz(nz*np-np+4)
              !      coca2 = cinfoz(nz*np-np+5)
              !      coca3 = cinfoz(nz*np-np+6)
              !      cocb1 = cinfoz(nz*np-np+7)
              !      cocb2 = cinfoz(nz*np-np+8)
              !      cocb3 = cinfoz(nz*np-np+9)
              !      cocc1 = cinfoz(nz*np-np+10)
              !      cocc2 = cinfoz(nz*np-np+11)
              !      cocc3 = cinfoz(nz*np-np+12)
              !      kdm=nint(cinfoz(nz*np-np+22))
              !      r1 = x1*x1 + y1*y1 + z1*z1
              !      xt = x1/(r1*r1*r1)
              !      yt = y1/(r1*r1*r1)
              !      zt = z1/(r1*r1*r1)
              !      ! calculate the jump condition of u
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
              !      wt = (unout - unin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      uzu = (u(i,j,k)-u(i,j,k-1)-tt)/hz+a1*wt
              !      ui(i,j,k)=ui(i,j,k)-2.0*w(i,j,k)*uzu


              !      ! calculate the jump condition of u1
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,u1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
              !      wt = (unout - unin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      uzu = (u1(i,j,k)-u1(i,j,k-1)-tt)/hz+a1*wt
              !      ui(i,j,k)=ui(i,j,k)+w1(i,j,k)*uzu
              !     
              !      ! calculate the jump condition of v
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
              !      wt = (vnout - vnin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      vzv = (v(i,j,k)-v(i,j,k-1)-tt)/hz+a1*wt
              !      vi(i,j,k)=vi(i,j,k)-2.0*w(i,j,k)*vzv


              !      ! calculate the jump condition of v1
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,v1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
              !      wt = (vnout - vnin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      vzv = (v1(i,j,k)-v1(i,j,k-1)-tt)/hz+a1*wt
              !      vi(i,j,k)=vi(i,j,k)+w1(i,j,k)*vzv
              !      
              !     ! calculate the jump condition of w
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
              !      wt = (wnout - wnin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      wzw = (w(i,j,k)-w(i,j,k-1)-tt)/hz+a1*wt
              !      wi(i,j,k)=wi(i,j,k)-2.0*w(i,j,k)*wzw


              !      ! calculate the jump condition of w1
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,w1,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
              !      wt = (wnout - wnin)*abs(coca3)
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      wzw = (w1(i,j,k)-w1(i,j,k-1)-tt)/hz+a1*wt
              !      wi(i,j,k)=wi(i,j,k)+w1(i,j,k)*wzw
              !      !         wt = 3.0d0*cc*cc*(xt*coca1+yt*coca2+zt*coca3)*abs(coca2)
              !      call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnout)
              !      call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
              !            ,phi,p2,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
              !            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,pnin,kdm,kzone)
              !      wt = (pnout - pnin)*abs(coca3)
              !      !         tt = -cc*cc/2.0d0/r1-1.0d0/r1
              !      a1 = (z1-z(k-1))/hz
              !      b1 = 1.d0 - a1
              !      pzp = (p2(i,j,k)-p2(i,j,k-1)-tt)/hz+a1*wt
              !      endif
              !      !         print *, 'pyp-1',i,j,a1,wt,pyp

                  else if(phi(i,j,k)*phi(i,j,k+1) <= 0.0 ) then

                     !         if(phi(i,j,k-1) .gt. 0.0) then
                     !           ui(i,j,k)= ui(i,j,k)-w(i,j,k)*(u(i,j,k)-u(i,j,k-1))/hz
                     !           vi(i,j,k)= vi(i,j,k)-w(i,j,k)*(v(i,j,k)-v(i,j,k-1))/hz
                     !           wi(i,j,k)= wi(i,j,k)-w(i,j,k)*(w(i,j,k)-w(i,j,k-1))/hz
                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(u(i,j,k)-u(i,j,k-1))/hz &
                           +ratio_t*w1(i,j,k)*(u1(i,j,k)-u1(i,j,k-1)) /hz
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(v(i,j,k)-v(i,j,k-1))/hz &
                           +ratio_t*w1(i,j,k)*(v1(i,j,k)-v1(i,j,k-1)) /hz
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(w(i,j,k)-w(i,j,k-1))/hz &
                           +ratio_t*w1(i,j,k)*(w1(i,j,k)-w1(i,j,k-1)) /hz
                     !         end if
                     pzp= (p2(i,j,k)-p2(i,j,k-1))/hz
                     !         print *, 'pyp+1',i,j,a1,wt,pyp

                     !       else if(phi(i,j)*phi(i-1,j) .le. 0.0 .or.
                     !    1          phi(i,j)*phi(i+1,j) .le. 0.0 ) then
                  else

                     !           ui(i,j,k)= ui(i,j,k)-w(i,j,k)*(u(i,j,k+1)-u(i,j,k-1))/hz2
                     !           vi(i,j,k)= vi(i,j,k)-w(i,j,k)*(v(i,j,k+1)-v(i,j,k-1))/hz2
                     !           wi(i,j,k)= wi(i,j,k)-w(i,j,k)*(w(i,j,k+1)-w(i,j,k-1))/hz2
                     ui(i,j,k)=ui(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(u(i,j,k+1)-u(i,j,k-1))/hz2 &
                           +ratio_t*w1(i,j,k)*(u1(i,j,k+1)-u1(i,j,k-1)) /hz2
                     vi(i,j,k)=vi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(v(i,j,k+1)-v(i,j,k-1))/hz2 &
                           +ratio_t*w1(i,j,k)*(v1(i,j,k+1)-v1(i,j,k-1)) /hz2
                     wi(i,j,k)=wi(i,j,k)-(1.0d0+ratio_t)*w(i,j,k)*(w(i,j,k+1)-w(i,j,k-1))/hz2 &
                           +ratio_t*w1(i,j,k)*(w1(i,j,k+1)-w1(i,j,k-1)) /hz2
                     pzp = (p2(i,j,k+1)-p2(i,j,k-1))/hz2

                  end if  ! (phi(i,j,k)*phi(i,j,k-1) <= 0.0)
                  !            else if(phi(i,j). gt. 0.0d0)  then             !# if 1

                  !       else
                  !            ui(i,j)= -(u(i,j)*(u(i+1,j)-u(i-1,j))/hx2
                  !    /          +v(i,j)*(u(i,j+1)-u(i,j-1))/hy2 )
                  !            vi(i,j)= -(u(i,j)*(v(i+1,j)-v(i-1,j))/hx2
                  !    /          +v(i,j)*(v(i,j+1)-v(i,j-1))/hy2 )
                  !      if (i==16.and.j==14.and.k==9) print *,'z',ui(i,j,k)
               end if  ! ( phi(i,j,k) > 0.d0 )
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)
               !external force
               if(phi(i,j,k) > 0.0) then
                  ui(i,j,k)=ui(i,j,k) + fe(t+0.5d0*dt,rnu,x(i),y(j),z(k))
                  vi(i,j,k)=vi(i,j,k) + ge(t+0.5d0*dt,rnu,x(i),y(j),z(k))
                  wi(i,j,k)=wi(i,j,k) + he(t+0.5d0*dt,rnu,x(i),y(j),z(k))

                  !         if(index2(i,j) .gt. 0) then

                  xp = x(i)
                  yp = y(j)
                  zp = z(k)
                  ne = 3
                  npp = 6
                  !           print *,'pxp',pxp,pyp
                  !           call inter_xya(m,n,npp,ne,hx,hy,t,dt,rnu,i,j,
                  !    1            phi,xp,yp,x,y,p,pxp,pyp,ww1,uw,vl)
                  !           print *,'pyp',pxp,pyp
                  !------------------------------------add pressure term
                  ui(i,j,k) = ui(i,j,k) -pxp
                  vi(i,j,k) = vi(i,j,k) -pyp
                  wi(i,j,k) = wi(i,j,k) -pzp
                  !         endif

                  !         if (i==28.and.j==32.and.k==32) print *,'one',is_out_ir
               end if
               !      if (i==16.and.j==14.and.k==9) print *,'p',ui(i,j,k)


            end do  !   i= 1,m-1
         end do  !   j= 1,n-1
      end do  !   k= 1,l-1
      !       stop
      !       write(*,*) 'ggg',ui(1,1,1)
      do k=0,l
         do j=0,n
            do i=0,m
               fj(i,j,k)=- ui(i,j,k)/rnu
               gj(i,j,k)= -vi(i,j,k)/rnu
               hj(i,j,k)= -wi(i,j,k)/rnu
            end do
         end do
      end do
      else         
      do k=0,l
         do j=0,n
            do i=0,m
               ui(i,j,k)=- fj(i,j,k)*rnu
               vi(i,j,k)= -gj(i,j,k)*rnu
               wi(i,j,k)= -hj(i,j,k)*rnu
            end do
         end do
      end do
      endif 
      ! here everything --CQ
      ! so inside there's no nonlinear term, no p gradient --CQ
      !----------------------Add the term taking the gradient of t
 call nsetimer_stop(2)
 call nsetimer_start(1)
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               !chnage the iteration of t
            ui(i,j,k)= -( ui(i,j,k) + ((1.0d0+ratio_t)**2*u(i,j,k)-ratio_t**2*u1(i,j,k))/((1.0d0+ratio_t)*dt))
            vi(i,j,k)= -( vi(i,j,k) + ((1.0d0+ratio_t)**2*v(i,j,k)-ratio_t**2*v1(i,j,k))/((1.0d0+ratio_t)*dt))
            wi(i,j,k)= -( wi(i,j,k) + ((1.0d0+ratio_t)**2*w(i,j,k)-ratio_t**2*w1(i,j,k))/((1.0d0+ratio_t)*dt))
               !      if (i==16.and.j==14.and.k==9) print *,'t',ui(i,j,k)
            end do
         end do
      end do
      ! u*/dt + ui = delta(u*)*rnu
      !       write(*,*) 'fff',ui(1,1,1)

      !------------------------------------# Add BC treatment:

      ! Neumann condition
      j = 0
      do k=1,l-1
         do i=1,m-1
            ui(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(u(i+1,j,k)-u(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(u1(i+1,j,k)-u1(i-1,j,k))/hx2 &
                  -(1.0+ratio_t)*w(i,j,k)*(u(i,j,k+1)-u(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(u1(i,j,k+1)-u1(i,j,k-1))/hz2
            vi(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(v(i+1,j,k)-v(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(v1(i+1,j,k)-v1(i-1,j,k))/hx2&
                  -(1.0+ratio_t)*w(i,j,k)*(v(i,j,k+1)-v(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(v1(i,j,k+1)-v1(i,j,k-1))/hz2
            wi(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(w(i+1,j,k)-w(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(w1(i+1,j,k)-w1(i-1,j,k))/hx2 &
                  -(1.0+ratio_t)*w(i,j,k)*(w(i,j,k+1)-w(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(w1(i,j,k+1)-w1(i,j,k-1))/hz2

            ui(i,j,k) = ui(i,j,k) - (p2(i+1,j,k)-p2(i-1,j,k))/hx2
            vi(i,j,k) = vi(i,j,k) - (p2(i,j+1,k)-p2(i,j,k))/hy
            wi(i,j,k) = wi(i,j,k) - (p2(i,j,k+1)-p2(i,j,k-1))/hz2

            ui(i,j,k)= -( ui(i,j,k) + ((1.0d0+ratio_t)**2*u(i,j,k)-ratio_t**2*u1(i,j,k))/((1.0d0+ratio_t)*dt))
            vi(i,j,k)= -( vi(i,j,k) + ((1.0d0+ratio_t)**2*v(i,j,k)-ratio_t**2*v1(i,j,k))/((1.0d0+ratio_t)*dt))
            wi(i,j,k)= -( wi(i,j,k) + ((1.0d0+ratio_t)**2*w(i,j,k)-ratio_t**2*w1(i,j,k))/((1.0d0+ratio_t)*dt))
         end do
      end do  !  k=1,l-1

      j = n
      do k=1,l-1
         do i=1,m-1
            ui(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(u(i+1,j,k)-u(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(u1(i+1,j,k)-u1(i-1,j,k))/hx2 &
                  -(1.0+ratio_t)*w(i,j,k)*(u(i,j,k+1)-u(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(u1(i,j,k+1)-u1(i,j,k-1))/hz2
            vi(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(v(i+1,j,k)-v(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(v1(i+1,j,k)-v1(i-1,j,k))/hx2&
                  -(1.0+ratio_t)*w(i,j,k)*(v(i,j,k+1)-v(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(v1(i,j,k+1)-v1(i,j,k-1))/hz2
            wi(i,j,k)= -(1.0+ratio_t)*u(i,j,k)*(w(i+1,j,k)-w(i-1,j,k))/hx2&
                  +ratio_t*u1(i,j,k)*(w1(i+1,j,k)-w1(i-1,j,k))/hx2 &
                  -(1.0+ratio_t)*w(i,j,k)*(w(i,j,k+1)-w(i,j,k-1))/hz2 &
                  +ratio_t*w1(i,j,k)*(w1(i,j,k+1)-w1(i,j,k-1))/hz2

            ui(i,j,k) = ui(i,j,k) - (p2(i+1,j,k)-p2(i-1,j,k))/hx2
            vi(i,j,k) = vi(i,j,k) - (p2(i,j,k)-p2(i,j-1,k))/hy
            wi(i,j,k) = wi(i,j,k) - (p2(i,j,k+1)-p2(i,j,k-1))/hz2

            ui(i,j,k)= -( ui(i,j,k) + ((1.0d0+ratio_t)**2*u(i,j,k)-ratio_t**2*u1(i,j,k))/((1.0d0+ratio_t)*dt))
            vi(i,j,k)= -( vi(i,j,k) + ((1.0d0+ratio_t)**2*v(i,j,k)-ratio_t**2*v1(i,j,k))/((1.0d0+ratio_t)*dt))
            wi(i,j,k)= -( wi(i,j,k) + ((1.0d0+ratio_t)**2*w(i,j,k)-ratio_t**2*w1(i,j,k))/((1.0d0+ratio_t)*dt))
         end do  !  i=1,m-1
      end do  !  k=1,l-1
      !       write(*,*) 'eee',ui(1,1,1)

call nsetimer_stop(1)
call nsetimer_start(2)
      !       print *,'uii',ui(31,57,49)
      !# --  Add the correction terms dimension by dimension -----------------
      !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
      !    1    ,vi(i,j,k),wi(i,j,k)

      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               !       if (i==24.and.j==26.and.k==24) cycle
               ! both inside and outside --CQ
               !# .................. between x_{i-1} and x_i ...........................

               ! apply correction terms to both inside and outside --CQ
               x1 = x(i) - hx
               x2 = x(i)
               x3 = x(i) + hx
               fl = phi(i-1,j,k)
               fc = phi(i,j,k)
               f3 = phi(i+1,j,k)
               
               !# -------------------------------- !# Irregular grid?
               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexx(i,j,k)
                  xi = cinfox(nx*np-np+1)
                  !find correction terms
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+18)=fuju
                  cinfox(nx*np-np+19)=fvjv
                  cinfox(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfox(nx*np-np+18)
                  fvjv=cinfox(nx*np-np+19)
                  fwjw=cinfox(nx*np-np+20)
                  endif
                  call xfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,xi,u,v,w,ui,vi,wi,u1,v1,w1,cinfox &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i-2,j,k) > 0.0) then
                     nx = indexx(i-1,j,k) + 1
                  else
                     nx = indexx(i-1,j,k)
                  end if
                  xi = cinfox(nx*np-np+1)
                  ilr = -1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+18)=fuju
                  cinfox(nx*np-np+19)=fvjv
                  cinfox(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfox(nx*np-np+18)
                  fvjv=cinfox(nx*np-np+19)
                  fwjw=cinfox(nx*np-np+20)
                  endif
                  call xfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,xi,u,v,w,ui,vi,wi,u1,v1,w1,cinfox &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                 !# End of (x_{i-1},x_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i-1,j,k) > 0.0) then
                     nx = indexx(i,j,k) + 1
                  else
                     nx = indexx(i,j,k)
                  end if
                  xi = cinfox(nx*np-np+1)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+18)=fuju
                  cinfox(nx*np-np+19)=fvjv
                  cinfox(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfox(nx*np-np+18)
                  fvjv=cinfox(nx*np-np+19)
                  fwjw=cinfox(nx*np-np+20)
                  endif
                  call xfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,xi,u,v,w,ui,vi,wi,u1,v1,w1,cinfox &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexx(i+1,j,k)
                  xi = cinfox(nx*np-np+1)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+18)=fuju
                  cinfox(nx*np-np+19)=fvjv
                  cinfox(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfox(nx*np-np+18)
                  fvjv=cinfox(nx*np-np+19)
                  fwjw=cinfox(nx*np-np+20)
                  endif
                  call xfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,xi,u,v,w,ui,vi,wi,u1,v1,w1,cinfox &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                   !# End of (x_i,x_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               !# ---------------------- between y_{j-1} and y_j ...........................

               x1 = y(j) - hy
               x2 = y(j)
               x3 = y(j) + hy
               fl = phi(i,j-1,k)
               fc = phi(i,j,k)
               f3 = phi(i,j+1,k)
               !# ------------     !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexy(i,j,k)
                  yj = cinfoy(nx*np-np+2)
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+18)=fuju
                  cinfoy(nx*np-np+19)=fvjv
                  cinfoy(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoy(nx*np-np+18)
                  fvjv=cinfoy(nx*np-np+19)
                  fwjw=cinfoy(nx*np-np+20)
                  endif

                  call yfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,yj,u,v,w,ui,vi,wi,u1,v1,w1,cinfoy &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)

               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i,j-2,k) > 0.0) then
                     nx = indexy(i,j-1,k) + 1
                  else
                     nx = indexy(i,j-1,k)
                  end if
                  yj = cinfoy(nx*np-np+2)
                  ilr = -1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+18)=fuju
                  cinfoy(nx*np-np+19)=fvjv
                  cinfoy(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoy(nx*np-np+18)
                  fvjv=cinfoy(nx*np-np+19)
                  fwjw=cinfoy(nx*np-np+20)
                  endif
                  call yfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,yj,u,v,w,ui,vi,wi,u1,v1,w1,cinfoy &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                 !# End of (x_{i-1},x_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i,j-1,k) > 0.0) then
                     nx = indexy(i,j,k) + 1
                  else
                     nx = indexy(i,j,k)
                  end if
                  yj = cinfoy(nx*np-np+2)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+18)=fuju
                  cinfoy(nx*np-np+19)=fvjv
                  cinfoy(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoy(nx*np-np+18)
                  fvjv=cinfoy(nx*np-np+19)
                  fwjw=cinfoy(nx*np-np+20)
                  endif
                  call yfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,yj,u,v,w,ui,vi,wi,u1,v1,w1,cinfoy &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexy(i,j+1,k)
                  yj = cinfoy(nx*np-np+2)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+18)=fuju
                  cinfoy(nx*np-np+19)=fvjv
                  cinfoy(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoy(nx*np-np+18)
                  fvjv=cinfoy(nx*np-np+19)
                  fwjw=cinfoy(nx*np-np+20)
                  endif
                  call yfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,yj,u,v,w,ui,vi,wi,u1,v1,w1,cinfoy &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                   !# End of (x_i,x_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !# ---------------------- between z_{j-1} and z_j ...........................

               x1 = z(k) - hz
               x2 = z(k)
               x3 = z(k) + hz
               fl = phi(i,j,k-1)
               fc = phi(i,j,k)
               f3 = phi(i,j,k+1)
               !# ------------     !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexz(i,j,k)
                  zk = cinfoz(nx*np-np+3)
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+18)=fuju
                  cinfoz(nx*np-np+19)=fvjv
                  cinfoz(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoz(nx*np-np+18)
                  fvjv=cinfoz(nx*np-np+19)
                  fwjw=cinfoz(nx*np-np+20)
                  endif

                  call zfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,zk,u,v,w,ui,vi,wi,u1,v1,w1,cinfoz &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)

               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i,j,k-2) > 0.0) then
                     nx = indexz(i,j,k-1) + 1
                  else
                     nx = indexz(i,j,k-1)
                  end if
                  zk = cinfoz(nx*np-np+3)
                  ilr = -1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+18)=fuju
                  cinfoz(nx*np-np+19)=fvjv
                  cinfoz(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoz(nx*np-np+18)
                  fvjv=cinfoz(nx*np-np+19)
                  fwjw=cinfoz(nx*np-np+20)
                  endif
                  call zfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,zk,u,v,w,ui,vi,wi,u1,v1,w1,cinfoz &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                 !# End of (z_{i-1},z_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i,j,k-1) > 0.0) then
                     nx = indexz(i,j,k) + 1
                  else
                     nx = indexz(i,j,k)
                  end if
                  zk = cinfoz(nx*np-np+3)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+18)=fuju
                  cinfoz(nx*np-np+19)=fvjv
                  cinfoz(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoz(nx*np-np+18)
                  fvjv=cinfoz(nx*np-np+19)
                  fwjw=cinfoz(nx*np-np+20)
                  endif
                  call zfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,zk,u,v,w,ui,vi,wi,u1,v1,w1,cinfoz &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexz(i,j,k+1)
                  zk = cinfoz(nx*np-np+3)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+18)=fuju
                  cinfoz(nx*np-np+19)=fvjv
                  cinfoz(nx*np-np+20)=fwjw
                  else           
                  fuju=cinfoz(nx*np-np+18)
                  fvjv=cinfoz(nx*np-np+19)
                  fwjw=cinfoz(nx*np-np+20)
                  endif
                  call zfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi &
                        ,x,y,z,zk,u,v,w,ui,vi,wi,u1,v1,w1,cinfoz &
                        ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)
               end if                   !# End of (z_i,z_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

           enddo
         enddo
       enddo
      !       write(*,*) 'bbb',ui(1,1,1)
call nsetimer_stop(2)
call nsetimer_start(3)
      if(ifzero .eqv. .true.) then
        return
       endif
      do k=0,l
         do j=0,n
            do i=0,m
               ui(i,j,k)= ui(i,j,k)/rnu
               vi(i,j,k)= vi(i,j,k)/rnu
               wi(i,j,k)= wi(i,j,k)/rnu
               !      if (i==16.and.j==14.and.k==9) print *,'rnu',ui(i,j,k)
            end do
         end do
      end do
      !       write(*,*) 'aaa',ui(1,1,1)
      ! u*/(dt*rnu) + ui = delta(u*)


         if(ifupdate .eqv. .true.) then
         solvernorm=0.0d0
         do k=1,l-1
            do j=1,n-1
               do i=1,m-1
               if(index2(i,j,k) .lt. 0) then
               solvernorm=solvernorm+abs(sqrt(ui(i,j,k)**2+&
           vi(i,j,k)**2+wi(i,j,k)**2)*h*h)/3
               endif
               enddo
             enddo
         enddo
           solvernorm=solvernorm/n2*h
          if(master) then
          write(6,*) "norm of velocity solver",solvernorm
          endif
          endif
      !#---- boundary condition
      ! Dirichlet condition
      t1 = t + dt

 
      i=0
      do k=0,l
         do j=0,n
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,xa,y(j),z(k))
           vi(i,j,k) = ve(t1,rnu,xa,y(j),z(k))
           wi(i,j,k) = we(t1,rnu,xa,y(j),z(k))
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do

      i=m
      do k=0,l
         do j=0,n
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,xb,y(j),z(k))
           vi(i,j,k) = ve(t1,rnu,xb,y(j),z(k))
           wi(i,j,k) = we(t1,rnu,xb,y(j),z(k))
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do

      k=0
      do j=0,n
         do i=0,m
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,x(i),y(j),za)
           vi(i,j,k) = ve(t1,rnu,x(i),y(j),za)
           wi(i,j,k) = we(t1,rnu,x(i),y(j),za)
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do

      k=l
      do j=0,n
         do i=0,m
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,x(i),y(j),zb)
           vi(i,j,k) = ve(t1,rnu,x(i),y(j),zb)
           wi(i,j,k) = we(t1,rnu,x(i),y(j),zb)
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do

      j=0
      do k=0,l
         do i=0,m
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,x(i),ya,z(k))
           vi(i,j,k) = ve(t1,rnu,x(i),ya,z(k))
           wi(i,j,k) = we(t1,rnu,x(i),ya,z(k))
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do

      j=n
      do k=0,l
         do i=0,m
           if(nse_debug .le. 8 .and. nse_debug .ne.0) then
           ui(i,j,k) = ue(t1,rnu,x(i),yb,z(k))
           vi(i,j,k) = ve(t1,rnu,x(i),yb,z(k))
           wi(i,j,k) = we(t1,rnu,x(i),yb,z(k))
           else
            ui(i,j,k) = 0.d0
            vi(i,j,k) = 0.d0
            wi(i,j,k) = 0.d0
           endif
         end do
      end do
      !----------------------------------solve u,v,w
      elmbda = -(1.0d0+2.0d0*ratio_t)/((1.0d0+ratio_t)*dt*rnu)
    !  mdimf = m + 1
    !  ndimf = n + 1
    !  mbdcnd = 1
    !  nbdcnd = 1
    !  lbdcnd = 1
      accept=1e-6
      dummy=0.0d0
      solvopt=1
      epsx=1.d0
      epsy=1.d0
      epsz=1.d0
     

        if(nse_debug .le.8 .and. nse_debug .ne. 0) then
        do j = 1, n-1; do i = 1, m-1
              ui(i,j,1 ) = ui(i,j,1 ) - ui(i,j,0)/(h*h)
              ui(i,j,l-1) = ui(i,j,l-1 )- ui(i,j,l)/(h*h)
        end do; end do
          
          
        do k = 1, l-1; do i = 1, m-1
              ui(i,1,k) = ui(i,1,k) - ui(i,0,k)/(h*h)
              ui(i,n-1,k) = ui(i,n-1,k ) - ui(i,n,k)/(h*h)
        end do; end do
     
     
        do k = 1, l-1; do j = 1, n-1
              ui(1,j,k) = ui(1,j,k ) - ui(0,j,k)/(h*h)
              ui(m-1,j,k) = ui(m-1,j,k ) - ui(m,j,k)/(h*h)
        end do; end do
        optx=0
        opty=0
        optz=0
        else
        optx=0
        opty=1
        optz=0
        endif
            call init_param(l-1,m-1,n-1,(l-1)*(m-1),(l-1)*(m-1)*(n-1),10000,dummy,accept,&
                 elmbda,1.0d0,h,1.9d0)

            call allocate_array(solvopt) 
            xso(:) = 0.d0
            call init_array(solvopt,epsx,epsy,epsz,-ui(1:m-1,1:n-1,1:l-1)*h*h,iv,xso,optx,opty,optz)
            call pb_iccg(ui(1:m-1,1:n-1,1:l-1),xso)  
      if(opty .eq. 1) then 
       do k=1,l-1
          do i=1,m-1
          ui(i,0,k)=ui(i,1,k) 
          ui(i,n,k)=ui(i,n-1,k)
          end do
       end do
       endif

        if(nse_debug .le.8 .and. nse_debug .ne. 0) then
        do j = 1, n-1; do i = 1, m-1
              vi(i,j,1 ) = vi(i,j,1 ) - vi(i,j,0)/(h*h)
              vi(i,j,l-1) = vi(i,j,l-1 )- vi(i,j,l)/(h*h)
        end do; end do
          
          
        do k = 1, l-1; do i = 1, m-1
              vi(i,1,k) = vi(i,1,k) - vi(i,0,k)/(h*h)
              vi(i,n-1,k) = vi(i,n-1,k ) - vi(i,n,k)/(h*h)
        end do; end do
     
     
        do k = 1, l-1; do j = 1, n-1
              vi(1,j,k) = vi(1,j,k ) - vi(0,j,k)/(h*h)
              vi(m-1,j,k) = vi(m-1,j,k ) - vi(m,j,k)/(h*h)
        end do; end do
        endif
            call init_param(l-1,m-1,n-1,(l-1)*(m-1),(l-1)*(m-1)*(n-1),10000,dummy,accept,&
                 elmbda,1.0d0,h,1.9d0)

            
            xso(:) = 0.d0
            call init_array(solvopt,epsx,epsy,epsz,-vi(1:m-1,1:n-1,1:l-1)*h*h,iv,xso,optx,opty,optz)
          call pb_iccg(vi(1:m-1,1:n-1,1:l-1),xso)

         if(opty .eq. 1) then 
          do k=1,l-1
          do i=1,m-1
          vi(i,0,k)=vi(i,1,k) 
          vi(i,n,k)=vi(i,n-1,k)
          end do
          end do
          endif

        if(nse_debug .le.8 .and. nse_debug .ne. 0) then
        do j = 1, n-1; do i = 1, m-1
              wi(i,j,1 ) = wi(i,j,1 ) - wi(i,j,0)/(h*h)
              wi(i,j,l-1) = wi(i,j,l-1 )- wi(i,j,l)/(h*h)
        end do; end do
          
          
        do k = 1, l-1; do i = 1, m-1
              wi(i,1,k) = wi(i,1,k) - wi(i,0,k)/(h*h)
              wi(i,n-1,k) = wi(i,n-1,k ) - wi(i,n,k)/(h*h)
        end do; end do
     
     
        do k = 1, l-1; do j = 1, n-1
              wi(1,j,k) = wi(1,j,k ) - wi(0,j,k)/(h*h)
              wi(m-1,j,k) = wi(m-1,j,k ) - wi(m,j,k)/(h*h)
        end do; end do
        endif
            call init_param(l-1,m-1,n-1,(l-1)*(m-1),(l-1)*(m-1)*(n-1),10000,dummy,accept,&
                 elmbda,1.0d0,h,1.9d0)

            
            xso(:) = 0.d0
            call init_array(solvopt,epsx,epsy,epsz,-wi(1:m-1,1:n-1,1:l-1)*h*h,iv,xso,optx,opty,optz)
          call pb_iccg(wi(1:m-1,1:n-1,1:l-1),xso)
          call deallocate_array(solvopt)
          if(opty .eq. 1) then 
          do k=1,l-1
          do i=1,m-1
          wi(i,0,k)=wi(i,1,k) 
          wi(i,n,k)=wi(i,n-1,k)
          end do
          end do
          endif
!   call hw3crt(xa,xb,m,mbdcnd,bdxa,bdxb,ya,yb,n,nbdcnd,bdya &
!         ,bdyb,za,zb,l,lbdcnd,bdza,bdzb,elmbda,mdimf &
!         ,ndimf,ui,pertrb,ierror,ww)
!     !       print *,'sui',pertrb,ierror
!     !       do k=1,l-1
!     !       do j=1,n-1
!     !       do i=1,m-1
!     !       write(32,'(3i6,2f12.6)') i,j,k,ui(i,j,k)
!     !    1    ,ue(0.d0,0.5d0,x(i),y(j),z(k))
!     !       end do
!     !       end do
!     !       end do
!     mbdcnd = 1
!     nbdcnd = 1
!     lbdcnd = 1

!     j=0
!     do k=0,l
!        do i=0,m
!           !          bdya(i,k) = uey(t1,rnu,x(i),ya,z(k))
!           bdya(i,k) = 0.d0
!        end do
!     end do

!     j=n
!     do k=0,l
!        do i=0,m
!           !          bdyb(i,k) = uey(t1,rnu,x(i),yb,z(k))
!           bdyb(i,k) = 0.d0
!        end do
!     end do
!        

!   call hw3crt(xa,xb,m,mbdcnd,bdxa,bdxb,ya,yb,n,nbdcnd,bdya &
!      ,bdyb,za,zb,l,lbdcnd,bdza,bdzb,elmbda,mdimf &
!      ,ndimf,vi,pertrb,ierror,ww)

!     mbdcnd = 1
!     nbdcnd = 1
!     lbdcnd = 1

!     j=0
!     do k=0,l
!        do i=0,m
!           !          bdya(i,k) = uey(t1,rnu,x(i),ya,z(k))
!           bdya(i,k) = 0.d0
!        end do
!     end do

!     j=n
!     do k=0,l
!        do i=0,m
!           !          bdyb(i,k) = uey(t1,rnu,x(i),yb,z(k))
!           bdyb(i,k) = 0.d0
!        end do
!     end do

!   call hw3crt(xa,xb,m,mbdcnd,bdxa,bdxb,ya,yb,n,nbdcnd,bdya &
!         ,bdyb,za,zb,l,lbdcnd,bdza,bdzb,elmbda,mdimf &
!         ,ndimf,wi,pertrb,ierror,ww)

!     !       call phiout(m,n,ui)
!     !---------------------------------check u,v,w

!   call checku(m,n,l,t1,x,y,z,phi,ui,e,i0,j0,k0,rnu)
!   call checkv(m,n,l,t1,x,y,z,phi,vi,e1,i1,j1,k1,rnu)
!   call checkw(m,n,l,t1,x,y,z,phi,wi,e2,i2,j2,k2,rnu)
  
!       write(*,'(a,3i6,f12.8,3i6,f12.8,3i6,f12.8,3i6)') 'err-1' &
!             ,m,n,l,e,i0,j0,k0,e1,i1,j1,k1,e2,i2,j2,k2
      
call nsetimer_stop(3)
      !       stop
      !#    intermediate velocity in incremental form
      !-----------------------------update u,v,w
      do k=0,l
         do j=0,n
            do i=0,m
               fw(i,j,k)=ui(i,j,k)
               gw(i,j,k)=vi(i,j,k)
               hw(i,j,k)=wi(i,j,k)
            end do
         end do
      end do

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bilinear_m here]
   subroutine bilinear_m(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
         ,x1,y1,z1,x,y,z,phi,u,v,w,uxj &
         ,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj,u1,v1,w1)
      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,phi(0:m,0:n,0:l)
      _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1
      _REAL_ uxj,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj,u1,v1,w1
      integer m,n,l,i,j,k
      _REAL_ htx,hty,htz,x0,xe,y0,ye,z0,ze

      i = nint((x1-xa)/hx)
      j = nint((y1-ya)/hy)
      k = nint((z1-za)/hz)

      htx = 2.0/hx
      hty = 2.0/hy
      htz = 2.0/hz
      x0 = 1.0 -( (x1-x(i))*htx -1.0)
      y0 = 1.0 -( (y1-y(j))*hty -1.0)
      z0 = 1.0 -( (z1-z(k))*htz -1.0)
      xe = 1.0 + ( (x1-x(i)) *htx -1.0)
      ye = 1.0 + ( (y1-y(j)) *hty -1.0)
      ze = 1.0 + ( (z1-z(k)) *htz -1.0)

      u1 = 0.125*(u(i,j,k)*x0*y0*z0+u(i+1,j,k)*xe*y0*z0 &
            +u(i,j+1,k)*x0*ye*z0+u(i+1,j+1,k)*xe*ye*z0 &
            +u(i,j,k+1)*x0*y0*ze+u(i+1,j,k+1)*xe*y0*ze &
            +u(i,j+1,k+1)*x0*ye*ze+u(i+1,j+1,k+1)*xe*ye*ze)
      v1 = 0.125*(v(i,j,k)*x0*y0*z0+v(i+1,j,k)*xe*y0*z0 &
            +v(i,j+1,k)*x0*ye*z0+v(i+1,j+1,k)*xe*ye*z0 &
            +v(i,j,k+1)*x0*y0*ze+v(i+1,j,k+1)*xe*y0*ze &
            +v(i,j+1,k+1)*x0*ye*ze+v(i+1,j+1,k+1)*xe*ye*ze)
      w1 = 0.125*(w(i,j,k)*x0*y0*z0+w(i+1,j,k)*xe*y0*z0 &
            +w(i,j+1,k)*x0*ye*z0+w(i+1,j+1,k)*xe*ye*z0 &
            +w(i,j,k+1)*x0*y0*ze+w(i+1,j,k+1)*xe*y0*ze &
            +w(i,j+1,k+1)*x0*ye*ze+w(i+1,j+1,k+1)*xe*ye*ze)

      ! only outside --CQ
      if(phi(i,j,k) > 0.0) then
         u1 = u1 + 0.125*x0*y0*z0 &
               *(uxj*(x(i)-x1)+uyj*(y(j)-y1)+uzj*(z(k)-z1))
         v1 = v1 + 0.125*x0*y0*z0 &
               *(vxj*(x(i)-x1)+vyj*(y(j)-y1)+vzj*(z(k)-z1))
         w1 = w1 + 0.125*x0*y0*z0 &
               *(wxj*(x(i)-x1)+wyj*(y(j)-y1)+wzj*(z(k)-z1))
      end if

      if(phi(i+1,j,k) > 0.0) then
         u1 = u1 + 0.125*xe*y0*z0 &
               *(uxj*(x(i+1)-x1)+uyj*(y(j)-y1)+uzj*(z(k)-z1))
         v1 = v1 + 0.125*xe*y0*z0 &
               *(vxj*(x(i+1)-x1)+vyj*(y(j)-y1)+vzj*(z(k)-z1))
         w1 = w1 + 0.125*xe*y0*z0 &
               *(wxj*(x(i+1)-x1)+wyj*(y(j)-y1)+wzj*(z(k)-z1))
      end if

      if(phi(i,j+1,k) > 0.0) then
         u1 = u1 + 0.125*x0*ye*z0 &
               *(uxj*(x(i)-x1)+uyj*(y(j+1)-y1)+uzj*(z(k)-z1))
         v1 = v1 + 0.125*x0*ye*z0 &
               *(vxj*(x(i)-x1)+vyj*(y(j+1)-y1)+vzj*(z(k)-z1))
         w1 = w1 + 0.125*x0*ye*z0 &
               *(wxj*(x(i)-x1)+wyj*(y(j+1)-y1)+wzj*(z(k)-z1))
      end if

      if(phi(i+1,j+1,k) > 0.0) then
         u1 = u1 + 0.125*xe*ye*z0 &
               *(uxj*(x(i+1)-x1)+uyj*(y(j+1)-y1)+uzj*(z(k)-z1))
         v1 = v1 + 0.125*xe*ye*z0 &
               *(vxj*(x(i+1)-x1)+vyj*(y(j+1)-y1)+vzj*(z(k)-z1))
         w1 = w1 + 0.125*xe*ye*z0 &
               *(wxj*(x(i+1)-x1)+wyj*(y(j+1)-y1)+wzj*(z(k)-z1))
      end if

      if(phi(i,j,k+1) > 0.0) then
         u1 = u1 + 0.125*x0*y0*ze &
               *(uxj*(x(i)-x1)+uyj*(y(j)-y1)+uzj*(z(k+1)-z1))
         v1 = v1 + 0.125*x0*y0*ze &
               *(vxj*(x(i)-x1)+vyj*(y(j)-y1)+vzj*(z(k+1)-z1))
         w1 = w1 + 0.125*x0*y0*ze &
               *(wxj*(x(i)-x1)+wyj*(y(j)-y1)+wzj*(z(k+1)-z1))
      end if

      if(phi(i+1,j,k+1) > 0.0) then
         u1 = u1 + 0.125*xe*y0*ze &
               *(uxj*(x(i+1)-x1)+uyj*(y(j)-y1)+uzj*(z(k+1)-z1))
         v1 = v1 + 0.125*xe*y0*ze &
               *(vxj*(x(i+1)-x1)+vyj*(y(j)-y1)+vzj*(z(k+1)-z1))
         w1 = w1 + 0.125*xe*y0*ze &
               *(wxj*(x(i+1)-x1)+wyj*(y(j)-y1)+wzj*(z(k+1)-z1))
      end if

      if(phi(i,j+1,k+1) > 0.0) then
         u1 = u1 + 0.125*x0*ye*ze &
               *(uxj*(x(i)-x1)+uyj*(y(j+1)-y1)+uzj*(z(k+1)-z1))
         v1 = v1 + 0.125*x0*ye*ze &
               *(vxj*(x(i)-x1)+vyj*(y(j+1)-y1)+vzj*(z(k+1)-z1))
         w1 = w1 + 0.125*x0*ye*ze &
               *(wxj*(x(i)-x1)+wyj*(y(j+1)-y1)+wzj*(z(k+1)-z1))
      end if

      if(phi(i+1,j+1,k+1) > 0.0) then
         u1 = u1 + 0.125*xe*ye*ze &
               *(uxj*(x(i+1)-x1)+uyj*(y(j+1)-y1)+uzj*(z(k+1)-z1))
         v1 = v1 + 0.125*xe*ye*ze &
               *(vxj*(x(i+1)-x1)+vyj*(y(j+1)-y1)+vzj*(z(k+1)-z1))
         w1 = w1 + 0.125*xe*ye*ze &
               *(wxj*(x(i+1)-x1)+wyj*(y(j+1)-y1)+wzj*(z(k+1)-z1))
      end if

      return
   end subroutine

   !*******************************************************************************
   



   !*******************************************************************************

   double precision function prod1(n,x,y)
   implicit none

   _REAL_ x(n), y(n)
   _REAL_ prod1
   integer i,n

   prod1 = 0.0d0
   do i=1,n
      prod1 = prod1 + x(i)*y(i)
   end do

   return
end function
!FIXME:                nse.f, line 2777: bad indentation level at end of subroutine copyvec


   !#---------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dextend here]
   subroutine dextend(m,n,p0)

      !#-- Extend u(0:m,0:n) to the boundary using 2-end extrapolation -----

      implicit none
      _REAL_ p0(0:m,0:n)
      integer m,n,i,j

      i=0
      do j=1,n-1
         p0(i,j)=3.*p0(i+1,j)-3.*p0(i+2,j)+p0(i+3,j)
      end do

      i=m
      do j=1,n-1
         p0(i,j)=3.*p0(i-1,j)-3.*p0(i-2,j)+p0(i-3,j)
      end do

      j=0
      do i=1,m-1
         p0(i,j)=3.*p0(i,j+1)-3.*p0(i,j+2)+p0(i,j+3)
      end do

      j=n
      do i=1,m-1
         p0(i,j)=3.*p0(i,j-1)-3.*p0(i,j-2)+p0(i,j-3)
      end do

      i=0
      j=0
      p0(i,j)=3.*p0(i+1,j+1)-3.*p0(i+2,j+2)+p0(i+3,j+3)

      i=0
      j=n
      p0(i,j)=3.*p0(i+1,j-1)-3.*p0(i+2,j-2)+p0(i+3,j-3)

      i=m
      j=0
      p0(i,j)=3.*p0(i-1,j+1)-3.*p0(i-2,j+2)+p0(i-3,j+3)

      i=m
      j=n
      p0(i,j)=3.*p0(i-1,j-1)-3.*p0(i-2,j-2)+p0(i-3,j-3)

      return
   end subroutine



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bound_info here]
   subroutine bound_info(m,n,l,np,n2,xa,xb,ya,yb,za,zb,hx,hy,hz &
         ,x,y,z,phi,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,kzone)

      !#-------------------------------------------------------------------------
      !--
      !--   This subroutine find the intersection of the boundary
      !--   represented by the grid level set function phi(i,j)
      !--   and the grid lines x=x(i) and y=y(j)
      !--
      !-- Output:
      !--   indexx(i,j) >=1, an inside irregular grid poits which means
      !--   that the boundary intersects with the grid line y=y(j).
      !--
      !--   cinfox: Let indexx(i,j) = n1, then cinfox(3*n1-2) is the
      !--   intersection; cinfox(3*n1-1)=cos, cinfox(3*n1)= sic,
      !--   where (coc,sic) is the unit normal derivative.
      !--
      !--   indexy: Has the similar meaning as indexx.
      !--   cinfoy: Has the similar meaning as cinfoy.
      !-------------------------------------------------------------------------
      !--
      implicit none
      
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l)&
              ,kzone(0:m,0:n,0:l) 
      _REAL_ cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2)
      integer m,n,l,np,n2
      _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,f1,f2,f3
      integer i,j,k,nx,ny,nz
      _REAL_ sibb,sibc,sicc,x1,x2,x3,xroot,yroot,zroot 

      do k=0,l
         do j=0,n
            do i=0,m
               indexx(i,j,k) = 0
               indexy(i,j,k) = 0
               indexz(i,j,k) = 0
            end do
         end do
      end do

      nx = 0
      ny = 0
      nz = 0

      do 10 k=1,l-1
         do 10 j=1,n-1
            do 10 i=1,m-1
               !.................. between x_{i-1} and x_i ...........................

               x1 = x(i) - hx
               x2 = x(i)
               x3 = x(i) + hx
               f1 = phi(i-1,j,k)
               f2 = phi(i,j,k)
               f3 = phi(i+1,j,k)

               !......                                         !# Irregular grid?
               ! interface --CQ
               if( f1 > 0.0 .and. f2 <= 0.0) then
                  nx = nx +1
                  indexx(i,j,k) = nx

                  xroot = root(x1,x2,x3,f1,f2,f3)
                  cinfox(nx*np-np+1) = xroot
                  cinfox(nx*np-np+2) = y(j)
                  cinfox(nx*np-np+3) = z(k)
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,xroot,y(j),z(k),x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfox(nx*np-np+4) = coca1
                  cinfox(nx*np-np+5) = coca2
                  cinfox(nx*np-np+6) = coca3
                  cinfox(nx*np-np+7) = cocb1
                  cinfox(nx*np-np+8) = cocb2
                  cinfox(nx*np-np+9) = cocb3
                  cinfox(nx*np-np+10) = cocc1
                  cinfox(nx*np-np+11) = cocc2
                  cinfox(nx*np-np+12) = cocc3
                  cinfox(nx*np-np+13) = sibb
                  cinfox(nx*np-np+14) = sicc
                  cinfox(nx*np-np+15) = sibc
                  cinfox(nx*np-np+16) = curv
                  cinfox(nx*np-np+22) = kzone(i,j,k)/100000
               
                  !              if (i==47.and.j==46.and.k==29) print *,'bound'
                  !    1           ,x1,x2,x3,xroot,nx
               end if

               !.................. between x_{i} and x_{i+1} ...........................

               if( f3 > 0.0 .and. f2 <= 0.0) then
                  nx = nx +1
                  if( indexx(i,j,k) == 0) then
                     indexx(i,j,k) = nx
                  end if

                  !-- Note: if indexx(i,j) > 0, then the boundary cuts between (x_{i-1},x_i)
                  !--       and (x_i,x_{i+1}). In this case, indexx(i,j) is the index of the
                  !--       intersection starting from (x_{i-1},x_i), and the increment of
                  !--       the index is 2, including both (x_{i-1},x_i) and (x_i,x_{i+1}).
                  !--
                  xroot = root(x2,x3,x1,f2,f3,f1)
                  cinfox(nx*np-np+1) = xroot
                  cinfox(nx*np-np+2) = y(j)
                  cinfox(nx*np-np+3) = z(k)
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,xroot,y(j),z(k),x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfox(nx*np-np+4) = coca1
                  cinfox(nx*np-np+5) = coca2
                  cinfox(nx*np-np+6) = coca3
                  cinfox(nx*np-np+7) = cocb1
                  cinfox(nx*np-np+8) = cocb2
                  cinfox(nx*np-np+9) = cocb3
                  cinfox(nx*np-np+10) = cocc1
                  cinfox(nx*np-np+11) = cocc2
                  cinfox(nx*np-np+12) = cocc3
                  cinfox(nx*np-np+13) = sibb
                  cinfox(nx*np-np+14) = sicc
                  cinfox(nx*np-np+15) = sibc
                  cinfox(nx*np-np+16) = curv
                  cinfox(nx*np-np+22) = kzone(i,j,k)/100000
               end if  ! ( f3 > 0.0 .and. f2 <= 0.0)

               !.................. between y_{j-1} and y_j ...........................

               x1 = y(j) - hy
               x2 = y(j)
               x3 = y(j) + hy
               f1 = phi(i,j-1,k)
               f2 = phi(i,j,k)
               f3 = phi(i,j+1,k)
               !......                                         !# Irregular grid?

               if( f1 > 0.0 .and. f2 <= 0.0) then
                  ny = ny +1
                  indexy(i,j,k) = ny
                  yroot = root(x1,x2,x3,f1,f2,f3)
                  cinfoy(ny*np-np+1) = x(i)
                  cinfoy(ny*np-np+2) = yroot
                  cinfoy(ny*np-np+3) = z(k)
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,x(i),yroot,z(k),x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfoy(ny*np-np+4) = coca1
                  cinfoy(ny*np-np+5) = coca2
                  cinfoy(ny*np-np+6) = coca3
                  cinfoy(ny*np-np+7) = cocb1
                  cinfoy(ny*np-np+8) = cocb2
                  cinfoy(ny*np-np+9) = cocb3
                  cinfoy(ny*np-np+10) = cocc1
                  cinfoy(ny*np-np+11) = cocc2
                  cinfoy(ny*np-np+12) = cocc3
                  cinfoy(ny*np-np+13) = sibb
                  cinfoy(ny*np-np+14) = sicc
                  cinfoy(ny*np-np+15) = sibc
                  cinfoy(ny*np-np+16) = curv
                  cinfoy(ny*np-np+22) = kzone(i,j,k)/100000
               end if

               !.................. between y_{j} and y_{j+1} ...........................

               if( f3 > 0.0 .and. f2 <= 0.0) then
                  ny = ny +1
                  if( indexy(i,j,k) == 0) then
                     indexy(i,j,k) = ny
                  end if

                  !-- Note: if indexy(i,j) > 0, then the boundary cuts between (y_{j-1},y_j)
                  !--       and (y_j,y_{j+1}). In this case, indexy(i,j) is the index of the
                  !--       intersection starting from (y_{j-1},y_j), and the increment of
                  !--       the index is 2, including both (y_{j-1},y_j) and (y_j,y_{j+1}).
                  !--
                  yroot = root(x2,x3,x1,f2,f3,f1)
                  cinfoy(ny*np-np+1) = x(i)
                  cinfoy(ny*np-np+2) = yroot
                  cinfoy(ny*np-np+3) = z(k)
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,x(i),yroot,z(k),x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfoy(ny*np-np+4) = coca1
                  cinfoy(ny*np-np+5) = coca2
                  cinfoy(ny*np-np+6) = coca3
                  cinfoy(ny*np-np+7) = cocb1
                  cinfoy(ny*np-np+8) = cocb2
                  cinfoy(ny*np-np+9) = cocb3
                  cinfoy(ny*np-np+10) = cocc1
                  cinfoy(ny*np-np+11) = cocc2
                  cinfoy(ny*np-np+12) = cocc3
                  cinfoy(ny*np-np+13) = sibb
                  cinfoy(ny*np-np+14) = sicc
                  cinfoy(ny*np-np+15) = sibc
                  cinfoy(ny*np-np+16) = curv
                  cinfoy(ny*np-np+22) = kzone(i,j,k)/100000
               end if  ! ( f3 > 0.0 .and. f2 <= 0.0)

               !.................. between z_{j-1} and z_j ...........................

               x1 = z(k) - hz
               x2 = z(k)
               x3 = z(k) + hz
               f1 = phi(i,j,k-1)
               f2 = phi(i,j,k)
               f3 = phi(i,j,k+1)
               !......                                         !# Irregular grid?

               if( f1 > 0.0 .and. f2 <= 0.0) then
                  nz = nz +1
                  indexz(i,j,k) = nz
                  zroot = root(x1,x2,x3,f1,f2,f3)
                  cinfoz(nz*np-np+1) = x(i)
                  cinfoz(nz*np-np+2) = y(j)
                  cinfoz(nz*np-np+3) = zroot
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,x(i),y(j),zroot,x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfoz(nz*np-np+4) = coca1
                  cinfoz(nz*np-np+5) = coca2
                  cinfoz(nz*np-np+6) = coca3
                  cinfoz(nz*np-np+7) = cocb1
                  cinfoz(nz*np-np+8) = cocb2
                  cinfoz(nz*np-np+9) = cocb3
                  cinfoz(nz*np-np+10) = cocc1
                  cinfoz(nz*np-np+11) = cocc2
                  cinfoz(nz*np-np+12) = cocc3
                  cinfoz(nz*np-np+13) = sibb
                  cinfoz(nz*np-np+14) = sicc
                  cinfoz(nz*np-np+15) = sibc
                  cinfoz(nz*np-np+16) = curv
                  cinfoz(nz*np-np+22) = kzone(i,j,k)/100000
               end if

               !.................. between z_{j} and z_{j+1} ...........................

               if( f3 > 0.0 .and. f2 <= 0.0) then
                  nz = nz +1
                  if( indexz(i,j,k) == 0) then
                     indexz(i,j,k) = nz
                  end if

                  zroot = root(x2,x3,x1,f2,f3,f1)
                  cinfoz(nz*np-np+1) = x(i)
                  cinfoz(nz*np-np+2) = y(j)
                  cinfoz(nz*np-np+3) = zroot
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,x(i),y(j),zroot,x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  cinfoz(nz*np-np+4) = coca1
                  cinfoz(nz*np-np+5) = coca2
                  cinfoz(nz*np-np+6) = coca3
                  cinfoz(nz*np-np+7) = cocb1
                  cinfoz(nz*np-np+8) = cocb2
                  cinfoz(nz*np-np+9) = cocb3
                  cinfoz(nz*np-np+10) = cocc1
                  cinfoz(nz*np-np+11) = cocc2
                  cinfoz(nz*np-np+12) = cocc3
                  cinfoz(nz*np-np+13) = sibb
                  cinfoz(nz*np-np+14) = sicc
                  cinfoz(nz*np-np+15) = sibc
                  cinfoz(nz*np-np+16) = curv
                  cinfoz(nz*np-np+22) = kzone(i,j,k)/100000
               end if
      10 continue

      return
   end subroutine

   !-----------------------------------------------------------------------

   double precision function root(x0,x1,x2,f0,f1,f2)
   implicit none
   
   _REAL_ x0,x1,x2,f0,f1,f2,b,c,a0,b0,c0,r1,r2,t

   !  function root returns the approximation of root between x0 and x1 if
   !  f0*f1 <=0 using quadratic interpolation.

   b = (f0-f1)/(x0-x1)
   c = f2 - f1 - b*(x2-x1)
   c = c/( (x2-x0)*(x2-x1))

   a0 = c
   b0 = b - c*(x0+x1)
   c0 = f1 -b*x1 + c*x0*x1

   if( a0 == 0) then
      root = -c0/b0
      return
   end if

   t = b0*b0-4.0*a0*c0

   !----- If t <=0, must be double root t is close to zero...............
   if(t <= 0.0 ) then
      root = -b0/(2.0*a0)
      return
   end if

   t = dsqrt(t)
   if(b0 >= 0.0) then
      r1 = (-b0-t)/(2.0*a0)
   else
      r1 = (-b0 + t)/(2.0*a0)
   end if

   r2 = c0/(a0*r1)

   if( x0 <= r1 + 1.0e-7 .and. r1 <= x1+1.0e-7) then
      root = r1
   else
      root = r2
   end if

   return
end function
!FIXME:                nse.f, line 3129: bad indentation level at end of subroutine bound_info



   !#***********************************************************************
   !#             SUBROUTINE CURVINF
   !#***********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine curvinf2 here]
   subroutine curvinf2(m,n,l,n2,np,index2,xa,xb,ya,yb,za,zb &
       ,h,hx,hy,hz,x,y,z,phi,cinfo,kzone)
      implicit none
 
     
      _REAL_ phi(0:m,0:n,0:l)
      integer index2(0:m,0:n,0:l)
      _REAL_ x(0:m),y(0:n) 
      _REAL_ z(0:l),cinfo(n2*np)
      integer        kzone(0:m,0:n,0:l)
       integer m,n,l,n2,np,i,i1,j,j1,k,k_1,k1,k2,k3,k4,k5,k6,k7,kk,nk,nn
       _REAL_ xa,xb,ya,yb,za,zb,h,hx,hy,hz,coca1,coca2,coca3,x1,y1,z1
       _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,sibb,sibc,sicc
       integer nn1      

      !------------- Index grid points ----------------------------------

      !       print *, '----------'
      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               nn = index2(i,j,k)
               nn1 = abs(nn)
               ! interface --CQ
               if(nn1 >= 1) then
                  call project(m-1,n-1,l-1,h,hx,hy,hz,i,j,k,x,y,z,phi,x1,y1,z1)
                  !               cinfo(nn1*(np-1) + 1) = x1
                  !               cinfo(nn1*(np-1)+ 2) = y1
                  cinfo( (nn1-1)*np + 1) = x1
                  cinfo( (nn1-1)*np + 2) = y1
                  cinfo( (nn1-1)*np + 3) = z1
                  call normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz &
                        ,x1,y1,z1,x,y,z,phi,coca1,coca2,coca3 &
                        ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3 &
                        ,sibb,sicc,sibc,curv)
                  !               cinfo(nn1*(np-1)+3) = coc
                  !               cinfo(nn1*(np-1)+4) = sic
                  !               cinfo(nn1*(np-1)+5) = curv
                  cinfo( (nn1-1)*np + 4) = coca1
                  cinfo( (nn1-1)*np + 5) = coca2
                  cinfo( (nn1-1)*np + 6) = coca3
                  cinfo( (nn1-1)*np + 7) = cocb1
                  cinfo( (nn1-1)*np + 8) = cocb2
                  cinfo( (nn1-1)*np + 9) = cocb3
                  cinfo( (nn1-1)*np +10) = cocc1
                  cinfo( (nn1-1)*np +11) = cocc2
                  cinfo( (nn1-1)*np +12) = cocc3
                  cinfo( (nn1-1)*np +13) = sibb
                  cinfo( (nn1-1)*np +14) = sicc
                  cinfo( (nn1-1)*np +15) = sibc
                  cinfo( (nn1-1)*np +16) = curv
                  !               write(57,'(3i6,3f12.6)') i,j,k,x(i),y(j),z(k)
                  !               write(57,'(3f12.6)') x1,y1,z1
                  !               write(57,'(3f12.6)') coc1,coc2,coc3
                  !               print *, 'interface', x(i),y(j)

                  if ( kzone(i,j,k) >= 100000 ) then
                     kk = kzone(i,j,k)/100000
                     cinfo((nn1-1)*np+17) = kk
                  else
                     i1 = int(sign(1.d0,x1-x(i)))
                     j1 = int(sign(1.d0,y1-y(j)))
                     k_1 = int(sign(1.d0,z1-z(k)))
                     k1 = kzone(i+i1,j+j1,k)/100000
                     k2 = kzone(i+i1,j,k)/100000
                     k3 = kzone(i,j+j1,k)/100000
                     k4 = kzone(i+i1,j+j1,k+k_1)/100000
                     k5 = kzone(i+i1,j,k+k_1)/100000
                     k6 = kzone(i,j+j1,k+k_1)/100000
                     k7 = kzone(i,j,k+k_1)/100000
                     nk = 0
                     ! check 2D code
                     if ( k1 >= 1 ) then
                        kk = k1
                        nk = 1
                     else if ( k2 >= 1 ) then
                        kk = k2
                        nk = 1
                     else if ( k3 >= 1 ) then
                        kk = k3
                        nk = 1
                     else if ( k4 >= 1 ) then
                        kk = k4
                        nk = 1
                     else if ( k5 >= 1 ) then
                        kk = k5
                        nk = 1
                     else if ( k6 >= 1 ) then
                        kk = k6
                        nk = 1
                     else if ( k7 >= 1 ) then
                        kk = k7
                        nk = 1
                     end if
                     if ( nk == 0 ) then
                        print *, "projection owner not found", &
                              i,j,i1,j1,k1,k2,k3,k4,k5,k6,k7
                        stop
                     end if
                     if ( nk > 1 ) then
                        print *, "projection owner undertermined", &
                              i,j,i1,j1,k1,k2,k3
                        stop
                     end if
                     !                   if ( match(kzone(i,j,k)) == 0 ) then
                     !                      print *, "domain not match",
                     !     1                      i,j,k,kzone(i,j,k)
                     !                      stop
                     !                   end if
                     cinfo((nn1-1)*np+17) = kk
                  end if  ! ( kzone(i,j,k) >= 100000 )
                  !               print *,'interface',nn1
                  !               print *, x1,y1,z1
                  !               print *, coc1,coc2,coc3
                  !               print *, curv
                  !               if (nn1==2) stop
               end if  ! (nn1 >= 1)
            end do  !   i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of ubrubroutine fjump here]
   subroutine fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,uj,vj,wj,x,y,z,cinfo &
         ,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,curv,x1,y1,z1,uji,ujd,ujdd &
         ,vji,vjd,vjdd,wji,wjd,wjdd,infodex,nx,kdm,phi)
     implicit none
     integer nl,ne
      parameter( nl = 16, ne=6)
     

      !# ********************************************************************
      !#   This routine evaluate the jump information [u], [u]', [u]'', at a
      !#   point (x1,y1) on the interface using an arc-length interpolation.
      !#
      !#   Input information:
      !#
      !#   Onput information:
      !#     uji,ujd,ujdd, are the jumps of [u], [u]', [u]''.
      !# -----------------------------------------------------------------------

      _REAL_ x(0:m),y(0:n),z(0:l),w3(nl),w3b(nl) &
            ,w1(ne,nl),cinfo(n2*np),uj(n2),vj(n2),wj(n2),w3c(nl) &
            ,sd(ne+1),ew(nl),uw(ne,ne),v(nl,nl),ew3(nl),w4(ne),ew2(nl) 
      integer index2(0:m,0:n,0:l)
      _REAL_  ss(nl),phi(0:m,0:n,0:l) 
      integer ix(nl),iy(nl),iz(nl)
      _REAL_ hx,hy,hz,xa,xb,ya,yb,za,zb,coca1,coca2,coca3      
      integer m,n,l,n2,np,i,j,k,i0,j0,k0,k1,kdx,nn,nn1,nsub,nn2,infodex
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,x1,y1,z1
      _REAL_ uji,ujd,ujdd,vji,vjd,vjdd,wji,wjd,wjdd,alf,coca,cocb,cocc
      _REAL_ ew4,h,pi,sy,sz,ts,w3d,x2,y2,z2,x3,y3,z3,dis
      integer i1,i2,if0,inf,j1,job,ii,i3,j3,k3,nx,kdm,kdm1

   !   ! common /radius/r0
      !       write(54,'(a,3i6)') 'fjump',i,j,k
      nn=nx+(infodex-1)*n2*2
      if(ifupdate .eqv. .true.) then
      h = min(hx,hy,hz)
      pi = 3.14159265358979d0
      alf=5.1d0*h
      if0 = int(alf/h) + 1

      i0 = int((x1-xa)/hx)
      j0 = int((y1-ya)/hy)
      k0 = int((z1-za)/hz)
      dis=4*h*h
      i1=i0
      j1=j0
      k1=k0
      do i3=0,1
         do j3=0,1
           do k3=0,1
           if((x1-x(i1+i3))**2+(y1-y(j1+j3))**2+(z1-z(k1+k3))**2 .le. dis .and.&
            phi(i1+i3,j1+j3,k1+k3)<=0) then
           dis=(x1-x(i1+i3))**2+(y1-y(j1+j3))**2+(z1-z(k1+k3))**2
           i0=i1+i3
           j0=j1+j3
           k0=k1+k3
           endif
           enddo
         enddo
      enddo

      w3 = 0.0d0
      w3b = 0.0d0
      w3c = 0.0d0
      w3d = 0.0d0
      sd = 0.0d0
      ew = 0.0d0
      uw = 0.0d0
      v = 0.0d0
      ew4 = 0.0d0
      ew3 = 0.0d0
      ew2 = 0.0d0
      w4 = 0.0d0
      
      !# ------------------------ Generate the matrix -----------------------

      nsub = 0
      do i2=0, if0       !# Starting from the center
         do i1 = i0-i2,i0+i2
            do j1 = j0-i2, j0+i2
               do k1 = k0-i2, k0+i2                 
                   if (i1<1 .or. j1<1 .or. k1<1) cycle
                  if(i1>m-1 .or. j1>n-1 .or. k1>l-1) cycle
                  nn1 = index2(i1,j1,k1)
                  !           if (i1==47.and.j1==46.and.k1==29) then
                  !             print *,'ff',nn1,uj(nn1),vj(nn1),wj(nn1)
                  !           end if
                  ! interface --CQ
                  ! was outside, now inside --CQ
                  if( (abs(i1-i0)+abs(j1-j0)+abs(k1-k0) == i2) &
                        .and. (nn1 < 0) ) then  !# outside irregular point
                     !            x2 = cinfo(nn1*(np-1) + 1)
                     !            y2 = cinfo(nn1*(np-1)+ 2)
                     !            coc2 = cinfo(nn1*(np-1)+3)
                     !            sic2 = cinfo(nn1*(np-1)+4)
                     nn1=abs(nn1)
                     kdm1=nint(cinfo(nn1*np-np+17))
                    if(kdm1 .eq. kdm .or. kdm1 .eq. 0) then 
              !     i3=0
              !     j3=0
              !     k3=0
              !    if(i1.ne.i0)  i3 = nint(sign(1.d0,x(i1)-x(i0)))
              !    if(j1.ne.j0)  j3 = nint(sign(1.d0,y(j1)-y(j0)))
              !    if(k1.ne.k0)  k3 = nint(sign(1.d0,z(k1)-z(k0)))
   !              if(phi(i1,j1,k1) >0) then
   !          if(phi(i1+1,j1,k1).le.0 .and. phi(i0+1,j0,k0) .gt.0) goto 500
   !          if(phi(i1,j1+1,k1).le.0 .and. phi(i0,j0+1,k0) .gt.0) goto 500    
   !          if(phi(i1,j1,k1+1).le.0 .and. phi(i0,j0,k0+1) .gt.0) goto 500
   !          if(phi(i1-1,j1,k1).le.0 .and. phi(i0-1,j0,k0) .gt.0) goto 500
   !          if(phi(i1,j1-1,k1).le.0 .and. phi(i0,j0-1,k0) .gt.0) goto 500    
   !          if(phi(i1,j1,k1-1).le.0 .and. phi(i0,j0,k0-1) .gt.0) goto 500
   !               else
   !          if(phi(i1-1,j1,k1).gt.0 .and. phi(i0+1,j0,k0) .gt.0) goto 500
   !          if(phi(i1,j1-1,k1).gt.0 .and. phi(i0,j0+1,k0) .gt.0) goto 500    
   !          if(phi(i1,j1,k1-1).gt.0 .and. phi(i0,j0,k0+1) .gt.0) goto 500
   !          if(phi(i1+1,j1,k1).gt.0 .and. phi(i0-1,j0,k0) .gt.0) goto 500
   !          if(phi(i1,j1+1,k1).gt.0 .and. phi(i0,j0-1,k0) .gt.0) goto 500    
   !          if(phi(i1,j1,k1+1).gt.0 .and. phi(i0,j0,k0-1) .gt.0) goto 500
   !      
   !               endif
                   x2 = cinfo( (nn1-1)*np + 1)
                     y2 = cinfo( (nn1-1)*np + 2)
                     z2 = cinfo( (nn1-1)*np + 3)
                     coca = cinfo( (nn1-1)*np + 4)
                     cocb = cinfo( (nn1-1)*np + 5)
                     cocc = cinfo( (nn1-1)*np + 6)

                     !            kdm = cinfo( (nn1-1)*np + 17)
                 !     do ii=1,nsub
                 !        i3=ix(ii)
                 !        j3=iy(ii)
                 !        k3=iz(ii)
                 !        nn2=abs(index2(i3,j3,k3))             
                 !    x3 = cinfo( (nn2-1)*np + 1)
                 !    y3 = cinfo( (nn2-1)*np + 2)
                  !   z3 = cinfo( (nn2-1)*np + 3)
                  !       call distan(x2,y2,z2,x3,y3,z3,dis)
                  !        if(dis < 0.1*h) then
                  !         goto 500
                  !        endif
                  !    enddo
                     !            if ( kdx == 0 ) then
                     !               kdx = kdm
                     !            else
                     !               if ( kdx /= kdm ) cycle
                     !            end if

                     call arc(x1,y1,z1,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
                           ,cocc1,cocc2,cocc3,x2,y2,z2,coca,cocb,cocc,sy,sz)
                     !            ts = abs(sy)
                     ts = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

                     !            r0 = 0.8d0
                     !            s1 = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                     !            alpha = acos(1.d0-s1/r0**2/2.d0)
                     !            s1 = r0*alpha
                     !            ts = abs(s1)
                     !            s = s1

                     !            write(54,'(6i6,3f12.6)') i0,j0,k0,i1,j1,k1,s,s1
                     !    1         ,phi(i1,j1,k1)
                     !            write(54,'(6f12.6)') x1,y1,z1,coca1,coca2,coca3
                     !            write(54,'(6f12.6)') x2,y2,z2,coca,cocb,cocc

                     !            write(*,'(a,3i6,3f12.6)') 'sss',i,j,k,ts,alf,s1
                     !            write(*,'(a,3i6,3f12.6)') 'sss',i,j,k,ss,alpha,s1
                     !            write(*,'(a,3i6,2f12.6)') 'sss',i,j,k,sy,sz
                     if(ts < alf .and. nsub < nl) then
                        !              ts = 1.0d0+dcos(pi*ts/alf)
                        ts = 1.d0
                        nsub = nsub + 1
                        ss(nsub) = ts
                        ix(nsub) = i1
                        iy(nsub) = j1
                        iz(nsub) = k1
                        w1(1,nsub) = 1.0d0*ts
                        w1(2,nsub) = sy*ts
                        w1(3,nsub) = sz*ts
                        w1(4,nsub) = sy*sy*ts*0.5d0
                        w1(5,nsub) = sz*sz*ts*0.5d0
                        w1(6,nsub) = sy*sz*ts
                        fjumpinfo((nn-1)*npfjump+1,nsub)=float(i1)
                        fjumpinfo((nn-1)*npfjump+2,nsub)=float(j1) 
                        fjumpinfo((nn-1)*npfjump+3,nsub)=float(k1) 
                        fjumpinfo((nn-1)*npfjump+4,nsub)=ts 
                     end if
                 endif ! kdm .eq. kdm1 
             end if  ! ( (abs(i1-i0)+abs(j1-j0)+abs(k1-k0) == i2)
            500 continue
         enddo
        enddo
       enddo

     enddo

      !----------------- Call least square routine --------------------------
      numfjump(nn)=nsub
      if ( nsub == 0 ) then
         print *,'fjump',i,j,k,x1,y1,z1,nsub
         nsub = 10000
         stop
      end if

      job = 11
      call dsvdc(w1,ne,ne,nsub,sd,ew,uw,ne,v,nl,w4,job,inf)

      do i1=1,ne
         if(abs(sd(i1)) > 1.0d-13) then
            ew(i1)=uw(1,i1)/sd(i1)
            ew2(i1)=uw(2,i1)/sd(i1)
            ew3(i1)=uw(3,i1)/sd(i1)
         else
            ew(i1)= 0.0d0
            ew2(i1)= 0.0d0
            ew3(i1)=0.0d0
         end if
      end do

      !---------------------- w3(i) is the solution -------------

      do i1=1,nsub
         w3(i1) = 0.0
         w3b(i1) = 0.0
         w3c(i1) = 0.0
         do j1=1,ne
            w3(i1) = w3(i1) + v(i1,j1)*ew(j1)
            w3b(i1) = w3b(i1) + v(i1,j1)*ew2(j1)
            w3c(i1) = w3c(i1) + v(i1,j1)*ew3(j1)
         end do
       fjumpinfo((nn-1)*npfjump+5,i1)=w3(i1) 
       fjumpinfo((nn-1)*npfjump+6,i1)=w3b(i1)
       fjumpinfo((nn-1)*npfjump+7,i1)=w3c(i1) 
      end do

      else 
       nsub=numfjump(nn)
      do i1=1,nsub
       ix(i1)= nint(fjumpinfo((nn-1)*npfjump+1,i1)) 
       iy(i1)= nint(fjumpinfo((nn-1)*npfjump+2,i1))
       iz(i1)= nint(fjumpinfo((nn-1)*npfjump+3,i1)) 
       ss(i1)= fjumpinfo((nn-1)*npfjump+4,i1) 
       w3(i1) =fjumpinfo((nn-1)*npfjump+5,i1) 
       w3b(i1)=fjumpinfo((nn-1)*npfjump+6,i1)
       w3c(i1)=fjumpinfo((nn-1)*npfjump+7,i1)
      end do
      endif
      !# ------------------ Summarize to get u_j' and u_j'' ------------------

      uji = 0.0d0
      ujd = 0.0d0
      ujdd = 0.0d0
      vji = 0.0d0
      vjd = 0.0d0
      vjdd = 0.0d0
      wji = 0.0d0
      wjd = 0.0d0
      wjdd = 0.0d0

      do  i2=1, nsub     !# Starting from the center
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)
         nn1 = abs(index2(i1,j1,k1))
       
         uji = uji + w3(i2)*uj(nn1)*ts
         ujd = ujd + w3b(i2)*uj(nn1)*ts
         ujdd = ujdd + w3c(i2)*uj(nn1)*ts
         vji = vji + w3(i2)*vj(nn1)*ts
         vjd = vjd + w3b(i2)*vj(nn1)*ts
         vjdd = vjdd + w3c(i2)*vj(nn1)*ts
         wji = wji + w3(i2)*wj(nn1)*ts
         wjd = wjd + w3b(i2)*wj(nn1)*ts
         wjdd = wjdd + w3c(i2)*wj(nn1)*ts
      enddo

      return
   endsubroutine



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


          subroutine gmres(m,n,l,mm,np,n2,n3,istart,indexx,indexy&
         ,indexz,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p&
         ,fw,gw,hw,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz&
         ,index2,cinfo,imax,x0,bf,tol,iter,error,kzone,nzone,ubc2,aret&
         ,fj,gj,hj)

          
          implicit none


          _REAL_  x(0:m),y(0:n),z(0:l),u(0:m,0:n,0:l),v(0:m,0:n,0:l)&
           ,w(0:m,0:n,0:l),p(0:m,0:n,0:l),p1(0:m,0:n,0:l)&
           ,fw(0:m,0:n,0:l),gw(0:m,0:n,0:l),hw(0:m,0:n,0:l)&
           ,fj(0:m,0:n,0:l),gj(0:m,0:n,0:l),hj(0:m,0:n,0:l)&
           ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l)&
           ,phi(0:m,0:n,0:l)
           integer  index2(0:m,0:n,0:l)
           _REAL_  cinfo(np*n2)&
           ,cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2)
           integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l)
           _REAL_ x0(n3),bf(n3),r(n3),vg(n3,mm),hg(n3,mm),s(n3)&
           ,x11(n3),vk(n3),ww(n3),hhj(mm,2),yy(n3)
           integer kzone(0:m,0:n,0:l)
           integer iter,m,n,l,np,n2,n3,istart,imax,mm,nzone
           integer i,i1,infon,j,j1,k,mi,i2
           _REAL_ serr,ser0,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,rnu,tol,error
           _REAL_ t,dr1,dtt,hit,ser1,st,aret(100),ubc2(n2),ser2

          iter = 0;hhj=0.d0

          do  j=1,imax

          ser0 = 0.0
          vg=0.0d0
          hg=0.0d0 
          yy=0.0d0
          ww=0.0d0

             call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
              ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
              ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz&
              ,index2,cinfo&
              ,x0,yy,bf,kzone,nzone,ubc2,aret,fj,gj,hj)


             call resid(n3,yy,bf,r)

          ser2 = dsqrt(prod1(n3,r,r))/dsqrt(prod1(n3,bf,bf))
       if(master) then
         write(6,*) ser2,tol 
       endif
             dr1 = dsqrt(prod1(n3,r,r))

    !         write(*,*)'in gmres, dr1=',dr1

             if(dr1 .le. 1e-18)return

             do i=1,n3
              vg(i,1) = r(i)/dr1
                s(i) = 0.0d0
             enddo

             s(1) = dr1

!   # ---------------- For i=1,2, ..., mm ------------------------------
             do i=1,mm-1

               iter = iter + 1

               do i1=1,n3
                   vk(i1) = vg(i1,i)
               enddo

             call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
              ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
              ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz&
              ,index2,cinfo&
              ,vk,ww,bf,kzone,nzone,ubc2,aret,fj,gj,hj)


   
!   # ---------------------------------------------------------------

          call resid(n3,ww,bf,r)
          ser2 = dsqrt(prod1(n3,r,r))/dsqrt(prod1(n3,bf,bf))
         if(master) then
          write(6,*) ser2,tol,iter 
         endif
          do k=1,i
             do i1=1,n3
                vk(i1) = vg(i1,k)
             enddo

             hg(k,i) = prod1(n3,ww,vk)

             do i1=1,n3
                ww(i1) = ww(i1) - hg(k,i)*vk(i1)
             enddo
          enddo

          hg(i+1,i) = dsqrt(prod1(n3,ww,ww))

          if( abs(hg(i+1,i)) .le. 1e-18) then
             infon = 1
          write(6,*), "infro=1,Got compleste basis of space"
          else
             infon = 0
             do i1=1,n3
                vg(i1,i+1) = ww(i1)/hg(i+1,i)
             enddo
          endif

 ! # ----------------------------------------------------------------

!  # ----- Apply J_1, j_2, ..., J_{i-1} on (hg_{1,i}, ..., hg_{i+1,i} ---
!  #   Suppose J_i =
!  #             | I 0 |    p = | cos(\alf)  sin(\alf) |
!  #               | 0 P |           | -sin(\alf) cos(\alf) |
!  #
!  #      cos(\alf) = hg(i,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)
!  #      sin(\alf) = hg(i+1,i)/sqrt(hg(i,i)^2 + hg(i+1,i)^2)

!  #------- Form J_i so that the (i+1)th component of J_i*h(:,i) is zero.

          do k=1,i-1
            hit =hhj(k,1)* hg(k,i) + hhj(k,2)*hg(k+1,i)
            hg(k+1,i) = -hhj(k,2)*hg(k,i) +hhj(k,1)*hg(k+1,i)
            hg(k,i) = hit
          enddo

          if(infon == 0) then

            dtt = dsqrt(hg(i,i)*hg(i,i)+hg(i+1,i)*hg(i+1,i))
            hhj(i,1) = hg(i,i)/dtt
            hhj(i,2) = hg(i+1,i)/dtt
            st =hhj(i,1)*s(i) + hhj(i,2)*s(i+1)
            s(i+1) = -hhj(i,2)*s(i) +hhj(i,1)*s(i+1)
            s(i) = st

            hit = hhj(i,1)* hg(i,i) +hhj(i,2)*hg(i+1,i)
            hg(i+1,i) = -hhj(i,2)*hg(i,i) + hhj(i,1)*hg(i+1,i)
            hg(i,i) = hit

            endif

          ser1 = abs(s(i+1))
      !    write(*,*)'***', i, ser1,infon


     !    if (ser1 .gt. 1.0d-15) then
     !       serr = abs((ser1-ser0)/ser1)
     !    else
     !       serr = 0.0
     !    endif

     !    if(ser1 .le. tol .or. serr .le. tol) then
     !       serr = tol - 1.0e-15
     !       mi = i
     !       goto 100
     !    endif


      !    enddo
             mi = i

!   #-------------- end loop for i=1,2,...mm-------------------------------

!   #--------------- Update(x_i,i) ----------------------------------------

       yy(mi) = s(mi)/(hg(mi,mi))

          do k=mi-1,1,-1
             yy(k) = s(k)
             do j1 = k+1,mi
               yy(k) = yy(k) - hg(k,j1)*yy(j1)
             enddo
             yy(k) = yy(k)/hg(k,k)
          enddo

          do i2=1,n3
            x11(i2) = x0(i2)
            do k=1,mi
               x11(i2) = x11(i2) + yy(k)*vg(i2,k)
            enddo
          enddo

          call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
           ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
           ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz&
           ,index2,cinfo&
           ,x11,ww,bf,kzone,nzone,ubc2,aret,fj,gj,hj)

          call resid(n3,ww,bf,r)
          ser2 = dsqrt(prod1(n3,r,r))/dsqrt(prod1(n3,bf,bf))
          
         if(master) then 
         write(6,*) ser2,tol,iter
         endif 
    !      error = s(mi+1)

    !      if( abs(s(mi+1)) .lt. tol .or. serr .lt. tol ) return
          if( abs(ser2) .lt. tol ) then
          mi=i
          goto 100
          endif

             ser0 = ser1
          enddo 
          mi=mm-1
    
      100 yy(mi) = s(mi)/(hg(mi,mi)+1.0d-14)
      !100    yy(mi) = s(mi)/hg(mi,mi)

      do k=mi-1,1,-1
         yy(k) = s(k)
         do j1 = k+1,mi
            yy(k) = yy(k) - hg(k,j1)*yy(j1)
         end do
         yy(k) = yy(k)/hg(k,k)
      end do

      !do i=1,nbnd
      do i=1,n3
         x11(i) = x0(i)
         do k=1,mi
            x11(i) = x11(i) + yy(k)*vg(i,k)
         end do
      end do
 
  !   write(*,*) 'XP: Following are used to test AVk=Vk+1Hk'
  !  
  !   write(*,*) ' XP: Following are tests for solution:'
 
  !stest=0.0d0     
  !   do i2=1,mi
  !    stest(i2)=s(i2)              
  !     do j2=1,mi
  !    stest(i2)=stest(i2)-hg(i2,j2)*yy(j2)              
  !     end do
  !  write(*,*) 'stest',i2,stest(i2)
  !   end do 
 
          call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
           ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
           ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz&
           ,index2,cinfo&
           ,x11,ww,bf,kzone,nzone,ubc2,aret,fj,gj,hj)
    !explicit residual

    call resid(n3,ww,bf,r)
    !  s(mi+1) = dsqrt(dot_product(r(1:ctn),r(1:ctn)))
    ! LiXiao changes Gmres Using relative error
    s(mi+1) =dsqrt(dot_product(r(1:n3),r(1:n3)))/dsqrt(dot_product(bf(1:n3),bf(1:n3)))
    !update x0 for potential restart.
     !-----Judge if residual OK reture, or back to loop 2 j=1 imax------------
 
    !  write(6,*) 'abs(s)=',abs(s(mi+1)),'tol=',tol,'j',j
    !  write(6,*) "relative:","RSD",RSD(ctn,r)/RSD0,"MAXD",MAXD(ctn,r)/MAXD0,"ABSD",ABSD(ctn,r)/ABSD0
    !  write(6,*) "absolute:","RSD",RSD(ctn,r),"MAXD",MAXD(ctn,r),"ABSD",ABSD(ctn,r)
    !  !!write(10,*) abs(s(mi+1)),iter

    do k=1,n3
       x0(k) = x11(k)
    end do
      if( abs(s(mi+1)) < tol ) then
      !if( RSD(ctn,r)/RSD0 < tol ) then
      !write(6,*) 'return'
      error=s(mi+1)
     return
     end if
     
   !avoid dead lock, when convergence can't be minimized, just restart.

   if (j > 4) then
   write(6,*) 'Possible dead lock in velocity, possible failure'
   ! write(6,*) 'j=',j,'res=',s(mi+1)
      return
     end if

      error = s(mi+1)

   end do  ! j=1,imax
   !write(6,*) 'j=',j,'imax=',imax
   return
             
          end subroutine

   !#*********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine resid here]
   subroutine resid(n,x,b,r)
      implicit none

      _REAL_ x(n), b(n), r(n)
     integer i,n
      do i=1,n
         r(i) = b(i) - x(i)
      end do

      return
   end subroutine

   !#*********************************************************************



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine iim_poisson_a here]
   subroutine iim_poisson_a(m,n,l,n_2,xa,xb,ya,yb,za,zb &
       ,phi,f,x,y,z,u,elmbda,ubc2,t,rnu)




      implicit none

      !#-------------------------------------------------------------------------
      !#
      !  Subroutine iim_poisson_a solvers an exterior Poisson/Helmholtz equation
      
      !       u_xx + u_yy + elmbda u = f(x,y)
      
      !  on a rectangular region with an arbitrary closed void region. The
      !  Dirichlet boundary condition is prescribed on the boundary of the void.
      
      
      !        -------------------------
      !        |                       |
      !        |        _-----\        |
      !        |       /       |       |
      !        |      /         \      |
      !        |      |          \     |
      !        |      |u=Given   |     |
      !        |      \          |     |
      !        |       \---------      | Given BC in the input format as
      !        | solution domain       | in Fish-pack
      !        |                       |
      !        -------------------------
      !#----------------------------------------------------------------
      !--  Input Parameters:
      !--
      !--    m, n:   The number of grid lines in x- and y- direction.

      !--           set it large enough (1240-- 12400), for example.
      !--    mm:    The number of previous vector used, (mm=80), for example.
      !--    a,b,c,d: The rectangular domain, a<=x<=b, c<=y<=d, which contains
      !--             the original irregular domain.
      !--    phi(0:m,0:n): The level set function. phi=0 is the boundary.
      !--                  phi <=0 is the original domain.
      !--    f(0:m,0:n): The right hand side of the Poisson equation when
      !--                phi<=0, f(i,j)=0 for phi >0.
      !--    x(0:m), y(0:n) The grid points in x- and y-directions.
      !--
      !--    
      !#--   ubc2(n2): the Dirichelt BC at the projections (augemented)
      !--
      !--  Parameters needed for Fish-pack routine: hwscrt:
      !--
      !--    mbdcnd,bda,bdb,nbdcnd,bdc,bdd,elmbda
      !--
      !--  Output Parameters:
      !--
      !--    u(0:m,0:n): The computed solution on the entire domain. The solution
      !--                to the original problem are at those points where
      !--                phi(i,j) >=0.
      !--    iter:       Number of iterations of GMRES iteration, which is also
      !--                the fast Poisson solver called.
      !--i
      !-- Working Parameters:
      !--
      !--    index(0:m,0:n), index2(0:m,0:n), work(2*(m+1)*(n+1))
      !--
      !-- Subroutine that users supplied:
      !--
      !--    function fji(x,y): The jump in the f(x,y), usually it is f(x,y)
      !--                       of the source on the boundary.
      !#-
      !--    function ubc(x,y): The Dirichlet boundary condition on the
      !--                       boundary.
      !--  The inputs of boundary condition of the rectangle is the same as
      !--  described in the fishpack subroutine: hwscrt.f
      !# ----------------------------------------------------------------------

      _REAL_ phi(0:m,0:n,0:l),f(0:m,0:n,0:l)&
            ,x(0:m),y(0:n),z(0:l)
      integer index(0:m,0:n,0:l)
      _REAL_  u(0:m,0:n,0:l)&
            ,ubc2(n_2)&
            ,phi_local(0:m+1,0:n+1,0:l+1),f_local(1:m-1,1:n-1,1:l-1) 
      integer iepsav(1:4,m*n*l),insas(0:m,0:n,0:l) 
      _REAL_ u_local(1:m-1,1:n-1,1:l-1)
      integer m,n,l,n_2
      _REAL_ xa,xb,ya,yb,za,zb,elmbda,t,rnu,accept,augctf
      _REAL_ augtol,epsin,epsout,h,idx,idy,idz,r
      integer i,iatmfirst,iatmlast,iaugtoltype,ibcopt,isolvopt
      integer j,k,n1,n2
      logical m_pbverbose
      _REAL_  norm,norm_i

      h = (xb -xa)/m
      epsin=1e0
      epsout=1e0
      accept=1e-4


      do k=0,l
         do j=0,n
            do i=0,m
               phi_local(i,j,k)=phi(i,j,k)
               if (phi(i,j,k) > 0) then
                  insas(i,j,k)= -1
               else
                  insas(i,j,k)= 1
               end if
            end do
         end do
      end do

      !#------------- Index grid points ----------------------------------


     n1 = 0
      n2 = 0
      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               index(i,j,k) = 0     !#  regular grid point
               if(phi(i,j,k) == 0.0) index(i,j,k) = 3
               if(phi(i,j,k) > 0.0) then
                  index(i,j,k)=5
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 4
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 4
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 4

               end if

               if(phi(i,j,k) < 0.0) then
                  index(i,j,k)=1
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 2
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 2
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 2

               end if

               if(index(i,j,k) == 4.or. index(i,j,k) == 3 .or. index(i,j,k)==2 & 
                ) then
                  n1 = n1 + 1
      
                  iepsav(1,n1) = i
                  iepsav(2,n1) = j
                  iepsav(3,n1) = k
      !            write(1,*) i,j,k,index(i,j,k),phi(i,j,k)
                  if(phi(i,j,k) > 0.0d0) then
                     n2 = n2 + 1

                  end if
               end if
      !         write(100,*) i,j,k,index(i,j,k),phi(i,j,k)
            end do  !   i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1
 !     stop     
 !#----------- Finish indexing grid points, below get interface information:
      !#--  projection: (x1,y1), cinfo(nn1,1), cinfo(nn1,2)
              
!--  normal direction: coc,sic: cinfo(nn1,3), cinfo(nn1,4)
      !--  curvature: curv, cinfo(nn1,5)
      !      stop
      !       write(*,*)'iim fore pressure n2='n, n2
      if (n_2/=n1) then
         write (*,*) "error in iim_poisson_a:number of &
               irregular points doesn't match" &
               ,n1,n_2
!         stop
      end if
      !        call curvinf1(m,n,l,n1,n2,index,xa,xb,ya,yb,za,zb
      !     1        ,h,h,h,x,y,z,phi,cinfo)
      iatmfirst = 0
      iatmlast = 0
      ibcopt = 1
      isolvopt = 1
      iaugtoltype = 0
      augctf = 0e0
      augtol = 1e-3
      m_pbverbose=.false.
      norm_i=0e0
      norm=0e0
      !          write(123,*) h

           
         

         

         ! k=0 and k=zm+1 faces
         if (nse_debug .le. 8 .and. nse_debug .ne.0) then
         do j = 1, n-1; do i = 1, m-1
            idx = abs(x(i)); idy = abs(y(j)); idz =abs(z(0))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(i,j,1 ) = f(i,j,1 ) - pe(rnu,t,x(i),y(j),z(0))/(h*h)
               u(i,j,0)=0 
            idz = abs(z(l))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(i,j,l-1) = f(i,j,l-1 )- pe(rnu,t,x(i),y(j),z(l))/(h*h)
               u(i,j,l)=0
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, l-1; do i = 1, m-1
            idx = abs(x(i)); idy  = abs(y(0)); idz  = abs(z(k))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(i,1,k) = f(i,1,k) - pe(rnu,t,x(i),y(0),z(k))/(h*h)
               u(i,0,k)=0
            idy = abs(y(n))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(i,n-1,k) = f(i,n-1,k ) - pe(rnu,t,x(i),y(n),z(k))/(h*h)
               u(i,n,k)=0
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, l-1; do j = 1, n-1
            idx = abs(x(0)); idy = abs(y(j)); idz = abs(z(k))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(1,j,k) = f(1,j,k ) - pe(rnu,t,x(0),y(j),z(k))/(h*h)
               u(0,j,k)=0
            idx = abs(x(m))
               r = sqrt(idx**2 + idy**2 + idz**2)
               f(m-1,j,k) = f(m-1,j,k ) - pe(rnu,t,x(m),y(j),z(k))/(h*h)
               u(m,j,k)=0 
         end do; end do
        endif

      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               u_local(i,j,k)= u(i,j,k)
               f_local(i,j,k)=f(i,j,k)

            end do
         end do
      end do
   
      
      call aug(x(0),y(0),z(0),m-1,n-1,l-1,phi,insas,n1 &
            ,iepsav,epsin,epsout,f_local &
            ,u_local,accept,h,iatmfirst &
            ,iatmlast,ibcopt &
            ,isolvopt,iaugtoltype,augctf,augtol,ubc2,m_pbverbose &
            ,norm_i,norm)
        if(nse_debug .le. 3 .and. nse_debug .ne.0) then
        return
        endif
      !            write(321,*) h
    
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               u(i,j,k)=u_local(i,j,k)

            end do
         end do
      end do
   end subroutine

   !#----------------------------------------------------------------




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexing2 here]
   subroutine indexing2(m,n,l,n2,index2,phi)
      implicit none

      _REAL_  phi(0:m,0:n,0:l)
      integer index2(0:m,0:n,0:l)
      integer m,n,l,n2,i,j,k,info

      !# ------------- Index grid points ----------------------------------

      n2 = 0
      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               index2(i,j,k) = 0
               !            if(phi(i,j) .eq. 0.0d0) index2(i,j) = 1

               !            if(phi(i,j) .gt. 0.0d0) then
               !            if(phi(i,j)*phi(i-1,j) .le. 0.0 .or.
               !    1                   phi(i,j)*phi(i+1,j) .le. 0.0) index2(i,j) = 1
               !            if(phi(i,j)*phi(i,j-1) .le. 0.0 .or.
               !    1                   phi(i,j)*phi(i,j+1) .le. 0.0) index2(i,j) = 1
               !            endif
     
       !     if(phi(i,j,k)>0) then          
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index2(i,j,k) = -1
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index2(i,j,k) = -1
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index2(i,j,k) = -1
      !          end if  
             info=1         
               if(phi(i,j,k)>0) then          
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .and. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) info=0
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .and. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) info=0
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .and. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) info=0
                end if           

               !            if(index2(i,j) .gt. 0) then
               !               n2 = n2 + 1
               !               index2(i,j) = n2
               !            endif
               if(index2(i,j,k) < 0 .and. info .eq. 1) then
                  n2 = n2 + 1
                  index2(i,j,k) = -n2
                     
               end if
            end do
         end do
      end do
      return
   end subroutine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexing2 here]
   subroutine mergepoint(m,n,l,phi,kzone)
      implicit none

      _REAL_  phi(0:m,0:n,0:l)
      integer kzone(0:m,0:n,0:l),index2(0:m,0:n,0:l)
      integer m,n,l,i,j,k,i0,j0,k0,i1,i2,j2,k2,if0,alf
       
     
     if0=1
     alf=1
      !# ------------- Index grid points ----------------------------------

    
      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               index2(i,j,k) = 0
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index2(i,j,k) = -1
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index2(i,j,k) = -1
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index2(i,j,k) = -1
                     
            end do
         end do
      end do
      do k=1,l-1
         do j=1,n-1
            do  i=1,m-1
               if(phi(i,j,k)>0 .and. index2(i,j,k)==-1) then
       do  i1=0,if0
         do i2=0,i1
            do j2=0,i1
               do k2=0,i1
                 if(i2**2+j2**2+k2**2<=alf) then
                 if(phi(i2+i,j2+j,k2+k) <=0  .and. phi(i-i2,j-j2,k-k2) <=0) then
                  kzone(i,j,k)=-1
                  endif
                 endif
             enddo
           enddo
         enddo
       enddo
     
               end if
            end do
         end do
      end do
   !  do k=1,l-1
   !     do j=1,n-1
   !        do  i=1,m-1
   !                if(kzone(i,j,k).eq.-1) then
   !                 phi(i,j,k)=-0.01
   !                 endif 
   !        end do
   !     end do
   !  end do
      return
   end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine interp3 here]
         subroutine interp3d(m,n,l,np,n2,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,t,rnu,phi,u,v,w,ft,fn,fm,x,y,z &
         ,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,curv,x1,y1,z1,uvout,vvout,wvout &
         ,utout,vtout,wtout,usout,vsout,wsout,uout,vout,wout,fuju,fvjv,fwjw,kzone)
      
      implicit none
    
      integer nl,ne
      parameter( nl=27, ne=10 )

      !#-------------------------------------------------------------------------
      !# Subroutine interp3 interpolates u^+ of u_{ij} at (x1,y1) which is a point
      !# on the interface. The method is the weighted least squares interpolation.
      !#-- m,n:          The number of grid in the x- and y-direction.
      !#-- x1,y1:  A point on the interface where the interpolation take place.
      !#-- u(0:m,0:n):   The grid function to be interpolated.
      !#-------------------------------------------------------------------------

      

      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l),u(0:m,0:n,0:l) &
            ,v(0:m,0:n,0:l),w(0:m,0:n,0:l)
       integer index2(0:m,0:n,0:l),kzone(0:m,0:n,0:l) 
       _REAL_ ft(n2),fn(n2),fm(n2),cinfo(np*n2),b(3,nl),w3(nl),w3a(nl),w3b(nl)&
            ,w3c(nl),w1(ne,nl),w2(nl,3),vl(nl,nl),sd(ne+1),ew1(nl),uw(ne,ne),w4(ne) &
            ,ew2(nl),ew3(nl),ss(nl),ew(nl)
       integer ix(nl),iy(nl),iz(nl)
       integer m,n,l,np,n2,i,j,k,infodex,nx
       _REAL_ hx,hy,hz,xa,xb,ya,yb,za,zb,t,rnu,coca1,coca2,coca3,cocb1,cocb2
       _REAL_ cocb3,cocc1,cocc2,cocc3,curv,x1,y1,z1,uvout,vvout,wvout
       _REAL_ usout,vsout,wsout,alf,ddfm,ddft,dfm,dfn,dft,fmv,fnv,dis       
       _REAL_ utout,vtout,wtout,ddfn,ftv,h,fuju,fvjv,fwjw,sibb,sicc,sibc
       integer i0,j0,k0,i1,i2,if0,ii,inf,j1,job,k1,nn,nsub,jj,kdm,i3,j3,k3
       _REAL_ ts,tt,ucort,unji,vcort,vnji,wcort,wnji,pi,uout,vout,wout

      !# --- Get jump conditions --------------------------------------------
      infodex=4
      nx=abs(index2(i,j,k))
     kdm=nint(cinfo((nx-1)*np+17)) 
      call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
            ,xa,xb,ya,yb,za,zb,ft,fn,fm,x,y,z,cinfo &
            ,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,curv,x1,y1,z1,ftv,dft,ddft &
            ,fnv,dfn,ddfn,fmv,dfm,ddfm,infodex,nx,kdm,phi)
      
      unji = ftv
      vnji = fnv
      wnji = fmv

      !       if (master) then
      !         r0 = 0.8d0
      !         r2 = (x1**2+y1**2+z1**2)
      !         r4 = r2*r2
      !         cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !         write(83,'(2f14.8)') ftv,-3.d0*cc*x1/(r4)
      !         write(83,'(2f14.8)') fnv,-3.d0*cc*y1/(r4)
      !         write(83,'(2f14.8)') fmv,-3.d0*cc*z1/(r4)
      !         print *,'fjump',x1,y1,z1,r2,-3.d0*cc*x1/r4,-3.d0*cc*y1/r4
      !    1      ,-3.d0*cc*z1/r4,ftv,fnv,fmv
      !       end if

      !       r0 = 0.8d0
      !       cc = (-2.0*rnu+sqrt(4.0*rnu*rnu-6.0*(r0-1.0/r0)+1.d-16))/3.d0

      !       print *,'interp3',unji,vnji,wnji
      !       print *,'interp3',-2.d0*cc*x1/r0**3,-2.d0*cc*y1/r0**3
      !    1    ,-2.d0*cc*z1/r0**3

      nn = abs(index2(i,j,k))
      unji = ft(nn)
      vnji = fn(nn)
      wnji = fm(nn)
     sibb=cinfo((nn-1)*np+13)
     sicc=cinfo((nn-1)*np+14)
     sibc=cinfo((nn-1)*np+15)
      if(ifupdate .eqv. .true.) then
      pi = 3.14159265358979d0
      h = min(hx,hy,hz)
      alf=5.1d0*h
      if0 = int(alf/h) + 1

      i0 = int((x1-xa)/h)
      j0 = int((y1-ya)/h)
      k0 = int((z1-za)/h)
      dis=4*h*h
      i1=i0
      j1=j0
      k1=k0
      do i3=0,1
         do j3=0,1
           do k3=0,1
           if((x1-x(i1+i3))**2+(y1-y(j1+j3))**2+(z1-z(k1+k3))**2 .le. dis .and.&
            phi(i1+i3,j1+j3,k1+k3)<=0) then
           dis=(x1-x(i1+i3))**2+(y1-y(j1+j3))**2+(z1-z(k1+k3))**2
           i0=i1+i3
           j0=j1+j3
           k0=k1+k3
           endif
           enddo
         enddo
      enddo



      !       unji = ftv/rnu
      !       vnji = fnv/rnu
      !       print *,'interp3',nn,ft(nn),fn(nn),fm(nn)

      !       uxj = unji*coca1
      !       uyj = unji*coca2
      !       uzj = unji*coca3
      !       vxj = vnji*coca1
      !       vyj = vnji*coca2
      !       vzj = vnji*coca3
      !       wxj = wnji*coca1
      !       wyj = wnji*coca2
      !       wzj = wnji*coca3

      !       uyyj = -dft/curv
      !       vyyj = -dfn/curv
      !       wyyj = -dfm/curv

      !       uyyjt = -unji*curv    !# = (ft * tau *curv )/rnu
      !       vyyjt = -vnji*curv
      !       wyyjt = -wnji*curv

      !       uxyjt = dft
      !       vxyjt = dfn
      !       wxyjt = dfm

      !       call bilinear_m(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz
      !    1    ,x1,y1,z1,x,y,z,phi,u,v,w,uxj
      !    1    ,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj,u1p,v1p,w1p)

      !       uxxjt = -uyyjt+(0.0*(u1p*coca1+v1p*coca2+w1p*coca3)*unji
      !    1    -fuju)/rnu
      !       vxxjt = -vyyjt+(0.0*(u1p*coca1+v1p*coca2+w1p*coca3)*vnji
      !    1    -fvjv)/rnu
      !       wxxjt = -wyyjt+(0.0*(u1p*coca1+v1p*coca2+w1p*coca3)*wnji
      !    1    -fwjw)/rnu

      ! uxz,uyz,uzz??

      !#------------------- Initialize ---------------------------------------

      do ii=1,nl
         do jj=1,nl
            vl(ii,jj) = 0.0d0
         end do
         do jj=1,ne
            w1(jj,ii) = 0.0d0
         end do
         w3a(ii) = 0.0
      end do
      w3a(2) = 1.0

      !#------------------------- Generate the matrix ------------------------

      nsub = 0
      do i2=0, if0       !# Starting from the center
         do i1 = i0-i2,i0+i2
            do j1 = j0-i2, j0+i2
               do  k1 = k0-i2, k0+i2
                   if (i1<1 .or. j1<1 .or. k1<1) cycle
                  if(i1>m-1 .or. j1>n-1 .or. k1>l-1) cycle
                  if( abs(i1-i0)+abs(j1-j0)+abs(k1-k0) == i2) then     !# if A
                     tt = (x(i1)-x1)*(x(i1)-x1)+(y(j1)-y1)* &
                           (y(j1)-y1)+(z(k1)-z1)*(z(k1)-z1)
                     ts = dsqrt(tt)
                  if(kdm .eq. kzone(i1,j1,k1)/100000 .or. phi(i1,j1,k1) .gt. 0.0) then !#if B 
         !         if( phi(i1,j1,k1) .le. 0.0) then !#if B 
              !     i3=0
              !     j3=0
              !     k3=0
              !    if(i1.ne.i0)  i3 = nint(sign(1.d0,x(i1)-x(i0)))
              !    if(j1.ne.j0)  j3 = nint(sign(1.d0,y(j1)-y(j0)))
              !    if(k1.ne.k0)  k3 = nint(sign(1.d0,z(k1)-z(k0)))
     !           if(phi(i1,j1,k1) >0) then
     !       if(phi(i1+1,j1,k1).le.0 .and. phi(i0+1,j0,k0) .gt.0) goto 600
     !       if(phi(i1,j1+1,k1).le.0 .and. phi(i0,j0+1,k0) .gt.0) goto 600    
     !       if(phi(i1,j1,k1+1).le.0 .and. phi(i0,j0,k0+1) .gt.0) goto 600
     !       if(phi(i1-1,j1,k1).le.0 .and. phi(i0-1,j0,k0) .gt.0) goto 600
     !       if(phi(i1,j1-1,k1).le.0 .and. phi(i0,j0-1,k0) .gt.0) goto 600    
     !       if(phi(i1,j1,k1-1).le.0 .and. phi(i0,j0,k0-1) .gt.0) goto 600
     !            else
     !       if(phi(i1-1,j1,k1).gt.0 .and. phi(i0+1,j0,k0) .gt.0) goto 600
     !       if(phi(i1,j1-1,k1).gt.0 .and. phi(i0,j0+1,k0) .gt.0) goto 600    
     !       if(phi(i1,j1,k1-1).gt.0 .and. phi(i0,j0,k0+1) .gt.0) goto 600
     !       if(phi(i1+1,j1,k1).gt.0 .and. phi(i0-1,j0,k0) .gt.0) goto 600
     !       if(phi(i1,j1+1,k1).gt.0 .and. phi(i0,j0-1,k0) .gt.0) goto 600    
     !       if(phi(i1,j1,k1+1).gt.0 .and. phi(i0,j0,k0-1) .gt.0) goto 600
     !   
     !             endif
                     if(ts < alf .and. nsub < nl &
                     ) then        !# if C
                        ts = 1.0d0+dcos(pi*ts/alf)
                       ts=1.d0 
                        nsub = nsub + 1
                        ix(nsub) = i1
                        iy(nsub) = j1
                        iz(nsub) = k1
                        ss(nsub) = ts
                        interpinfo((nn-1)*npinterp+1,nsub)=float(i1)
                        interpinfo((nn-1)*npinterp+2,nsub)=float(j1)
                        interpinfo((nn-1)*npinterp+3,nsub)=float(k1)
                        interpinfo((nn-1)*npinterp+4,nsub)=ts
                        call trans(x1,y1,z1,x(i1),y(j1),z(k1) &
                              ,w2(nsub,1),w2(nsub,2),w2(nsub,3),coca1,coca2,coca3 &
                              ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3)
                        !          w2(nsub,1) = x(i1)
                        !          w2(nsub,2) = y(j1)
                        !          w2(nsub,3) = z(k1)
                        w1(1,nsub) = ts
                        w1(2,nsub) = w2(nsub,1)*ts
                        w1(3,nsub) = w2(nsub,2)*ts
                        w1(4,nsub) = w2(nsub,3)*ts
                        w1(5,nsub) = 0.5*w2(nsub,1)*w2(nsub,1)*ts
                        w1(6,nsub) = 0.5*w2(nsub,2)*w2(nsub,2)*ts
                        w1(7,nsub) = 0.5*w2(nsub,3)*w2(nsub,3)*ts
                        w1(8,nsub) = w2(nsub,1)*w2(nsub,2)*ts
                        w1(9,nsub) = w2(nsub,1)*w2(nsub,3)*ts
                        w1(10,nsub) = w2(nsub,3)*w2(nsub,2)*ts
                        interpinfo((nn-1)*npinterp+5,nsub)=w2(nsub,1)
                        interpinfo((nn-1)*npinterp+6,nsub)=w2(nsub,2)
                        interpinfo((nn-1)*npinterp+7,nsub)=w2(nsub,3)
                        !          if ( nsub==nl ) print *,'i2',i2
                         
                        !---------- End coordinates transformation -------------------------------

                     end if                                        !# end if C
                     end if                                       !# end if B
                  end if                                        !# end if A.
               600 continue
            enddo
          enddo
        enddo
      enddo
      numinterp(nn)=nsub
      !#----------------- Call least square routine --------------------------

      job = 11
      call dsvdc(w1,ne,ne,nsub,sd,ew,uw,ne,vl,nl,w4,job,inf)

      do i1=1,ne
         if(abs(sd(i1)) > 1.0d-13) then
            ew(i1)=uw(1,i1)/sd(i1)
            ew1(i1)=uw(2,i1)/sd(i1)
            ew2(i1)=uw(3,i1)/sd(i1)
            ew3(i1)=uw(4,i1)/sd(i1)
         else
            ew(i1)= 0.0d0
            ew1(i1)=0.0d0
            ew2(i1)= 0.0d0
            ew3(i1)= 0.0d0
         end if
      end do

      do i1=1,nsub
         w3(i1) = 0.0
         w3a(i1) = 0.0
         w3b(i1) = 0.0
         w3c(i1) = 0.0
         do j1=1,ne
            w3(i1) = w3(i1) + vl(i1,j1)*ew(j1)
            w3a(i1) = w3a(i1) + vl(i1,j1)*ew1(j1)
            w3b(i1) = w3b(i1) + vl(i1,j1)*ew2(j1)
            w3c(i1) = w3c(i1) + vl(i1,j1)*ew3(j1)
         end do
       interpinfo((nn-1)*npinterp+8,i1)=w3(i1)
       interpinfo((nn-1)*npinterp+9,i1)=w3a(i1)
       interpinfo((nn-1)*npinterp+10,i1)=w3b(i1)
       interpinfo((nn-1)*npinterp+11,i1)=w3c(i1)
      end do
      else 
      nsub=numinterp(nn)
      do i1=1,nsub
       ix(i1)=nint(interpinfo((nn-1)*npinterp+1,i1))
       iy(i1)=nint(interpinfo((nn-1)*npinterp+2,i1))
       iz(i1)=nint(interpinfo((nn-1)*npinterp+3,i1))
       ss(i1)=interpinfo((nn-1)*npinterp+4,i1)
       w2(i1,1)=interpinfo((nn-1)*npinterp+5,i1)
       w2(i1,2)=interpinfo((nn-1)*npinterp+6,i1)
       w2(i1,3)=interpinfo((nn-1)*npinterp+7,i1)
       w3(i1)=interpinfo((nn-1)*npinterp+8,i1)
       w3a(i1)=interpinfo((nn-1)*npinterp+9,i1)
       w3b(i1)=interpinfo((nn-1)*npinterp+10,i1)
       w3c(i1)=interpinfo((nn-1)*npinterp+11,i1)
      enddo
      endif 
    
      !#-------- Begin to compute the correction term for u^{-} ------------

      uout=0.0
      vout=0.0
      wout=0.0
      uvout = 0.0
      vvout = 0.0
      wvout = 0.0
      utout = 0.0
      vtout = 0.0
      wtout = 0.0
      usout = 0.0
      vsout = 0.0
      wsout = 0.0
      ucort = 0.0
      vcort = 0.0
      wcort = 0.0

      !       uxout = 0.0
      !       uyout = 0.0
      !       uzout = 0.0
      !       vxout = 0.0
      !       vyout = 0.0
      !       vzout = 0.0
      !       wxout = 0.0
      !       wyout = 0.0
      !       wzout = 0.0



      do i2=1, nsub
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)

         if (phi(i1,j1,k1)>0.0d0) then
           b(1,i2)=u(i1,j1,k1)*ts
           b(2,i2)=v(i1,j1,k1)*ts
           b(3,i2)=w(i1,j1,k1)*ts
            else
           if(solverorder .eq. 1) then
           b(1,i2)=ts*(u(i1,j1,k1)+w2(i2,1)*unji)
           b(2,i2)=ts*(v(i1,j1,k1)+w2(i2,1)*vnji)
           b(3,i2)=ts*(w(i1,j1,k1)+w2(i2,1)*wnji)
           else
           b(1,i2)=ts*(u(i1,j1,k1)+w2(i2,1)*unji&
            +0.5*w2(i2,1)*w2(i2,1)*(unji*(sibb+sicc)+fuju)& 
            -0.5*w2(i2,2)*w2(i2,2)*unji*sibb-0.5*w2(i2,3)*w2(i2,3)*unji*sicc&   
            +dft*w2(i2,1)*w2(i2,2)+ddft*w2(i2,1)*w2(i2,3)-unji*sibc*w2(i2,2)*w2(i2,3))
            b(2,i2)=ts*(v(i1,j1,k1)+w2(i2,1)*vnji&
            +0.5*w2(i2,1)*w2(i2,1)*(vnji*(sibb+sicc)+fvjv)& 
            -0.5*w2(i2,2)*w2(i2,2)*vnji*sibb-0.5*w2(i2,3)*w2(i2,3)*vnji*sicc&   
            +dfn*w2(i2,1)*w2(i2,2)+ddfn*w2(i2,1)*w2(i2,3)-vnji*sibc*w2(i2,2)*w2(i2,3))
             b(3,i2)=ts*(w(i1,j1,k1)+w2(i2,1)*wnji&
             +0.5*w2(i2,1)*w2(i2,1)*(wnji*(sibb+sicc)+fwjw)& 
             -0.5*w2(i2,2)*w2(i2,2)*wnji*sibb-0.5*w2(i2,3)*w2(i2,3)*wnji*sicc&   
             +dfm*w2(i2,1)*w2(i2,2)+ddfm*w2(i2,1)*w2(i2,3)-wnji*sibc*w2(i2,2)*w2(i2,3))
             endif
           endif  
         uout=uout+w3(i2)*b(1,i2)
         vout=vout+w3(i2)*b(2,i2)
         wout=wout+w3(i2)*b(3,i2)
         uvout= uvout + w3a(i2)*b(1,i2)
         vvout= vvout + w3a(i2)*b(2,i2)
         wvout= wvout + w3a(i2)*b(3,i2)
         utout= utout + w3b(i2)*b(1,i2)
         vtout= vtout + w3b(i2)*b(2,i2)
         wtout= wtout + w3b(i2)*b(3,i2)
         usout= usout + w3c(i2)*b(1,i2)
         vsout= vsout + w3c(i2)*b(2,i2)
         wsout= wsout + w3c(i2)*b(3,i2)
         !                uxout= uxout + w3(i2)*u(i1,j1,k1)*ts
         !                uyout= uyout + w3b(i2)*u(i1,j1,k1)*ts
         !                uzout= uzout + w3c(i2)*u(i1,j1,k1)*ts
         !                vxout= vxout + w3(i2)*v(i1,j1,k1)*ts
         !                vyout= vyout + w3b(i2)*v(i1,j1,k1)*ts
         !                vzout= vzout + w3c(i2)*v(i1,j1,k1)*ts
         !                wxout= wxout + w3(i2)*w(i1,j1,k1)*ts
         !                wyout= wyout + w3b(i2)*w(i1,j1,k1)*ts
         !                wzout= wzout + w3c(i2)*w(i1,j1,k1)*ts

         !                print *,'vout',i1,j1,k1,u(i1,j1,k1),cc*x(i1)
         !    1             /(x(i1)**2+y(j1)**2+z(k1)**2)

         ! temporarily use first-order correction
    !    if( phi(i1,j1,k1) <= 0.0d0 ) then
    !       !                ucort = ucort - w3(i2)*ts*( w2(i2,1)*(unji+
    !       !    1                w2(i2,1)*uxxjt/2.0 + w2(i2,2)*uxyjt ) +
    !       !    2                w2(i2,2)*w2(i2,2)*uyyjt/2.0 )
    !       !                vcort = vcort - w3(i2)*ts*( w2(i2,1)*(vnji+
    !       !    1                w2(i2,1)*vxxjt/2.0 + w2(i2,2)*vxyjt ) +
    !       !    2                w2(i2,2)*w2(i2,2)*vyyjt/2.0 )
    !       ucort = ucort - w3(i2)*ts*w2(i2,1)*unji
    !       vcort = vcort - w3(i2)*ts*w2(i2,1)*vnji
    !       wcort = wcort - w3(i2)*ts*w2(i2,1)*wnji
    !    end if

      enddo

      !       uvout1 = uxout*coc1+uyout*coc2+uzout*coc3
      !       vvout1 = vxout*coc1+vyout*coc2+vzout*coc3
      !       wvout1 = wxout*coc1+wyout*coc2+wzout*coc3

      !       r0 = 0.8d0
      !       cc = (-2.0*rnu+sqrt(4.0*rnu*rnu-6.0*(r0-1.0/r0)+1.d-16))/3.d0
      !       print *, 'vout',cc*x1/r0**3,cc*y1/r0**3,cc*z1/r0**3
      !       print *, 'vout',ucort,vcort,wcort
      !       print *, 'vout',uvout1,vvout1,wvout1
  !   uvout = uvout - ucort
  !   vvout = vvout - vcort
  !   wvout = wvout - wcort
      !       print *, 'vout',uvout,vvout,wvout
      !       print *, 'tout',utout,vtout,wtout
      !       print *, 'sout',usout,vsout,wsout
      !       print *, 'vouta',-cc*x1/r0**3,-cc*y1/r0**3,-cc*z1/r0**3
      !       print *, 'touta',0.d0
      !    1    ,cc*z1/r0**2/dsqrt(y1**2+z1**2)
      !    1    ,-cc*y1/r0**2/dsqrt(y1**2+z1**2)
      !       print *, 'souta',-cc*dsqrt(y1**2+z1**2)/r0**3
      !    1    ,cc*x1*y1/r0**3/dsqrt(y1**2+z1**2)
      !    1    ,cc*x1*z1/r0**3/dsqrt(y1**2+z1**2)
      !       stop

      return
   end subroutine

   !#*********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trans here]
   subroutine trans(x,y,z,x1,y1,z1,x3,y3,z3,coca1,coca2,coca3 &
         ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3)
      implicit none
     _REAL_ x,y,z,x1,y1,z1,x3,y3,z3,coca1,coca2,coca3,cocb1,cocb2
     _REAL_ cocb3,cocc1,cocc2,cocc3

      !       call gettan(coc1,coc2,coc3,coca1,coca2,coca3
      !    1    ,cocb1,cocb2,cocb3)

      x3 = (x1-x)*coca1 + (y1-y)*coca2 + (z1-z)*coca3  !#(xi,eta) at(x2,y2)
      y3 = (x1-x)*cocb1 + (y1-y)*cocb2 + (z1-z)*cocb3
      z3 = (x1-x)*cocc1 + (y1-y)*cocc2 + (z1-z)*cocc3

      !       write(19,*) coc1,coc2,coc3,coca1,coca2,coca3
      !    1     ,coc1*coc1+coc2*coc2+coc3*coc3
      !    1     ,coc1*coca1+coc2*coca2+coc3*coca3
      !    1     ,coc1*cocb1+coc2*cocb2+coc3*cocb3
      !    1     ,cocb1*coca1+cocb2*coca2+cocb3*coca3

      return
   end subroutine




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine interp4 here]
   subroutine interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
         ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
      implicit none
      integer nl,ne
      parameter( nl = 16, ne = 10 )

      !--------------------------------------------------------------------------
      !-- Subroutine interp3 interpolates u_n^- and u_n^+ of the solution
      !-- of Poisson equation with jumps along interface at (x1,y1) which is a
      !-- point on the interface. The method is the weighted least squares
      !-- interpolation.
      !--
      !-- m,n:          The number of grid in the x- and y-direction.
      !-- x1,y1:        A point on the interface where the interpolation
      !--               take place.
      !-- u(0:m,0:n):   The grid function to be interpolated.
      !--------------------------------------------------------------------------

     


      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l),u(0:m,0:n,0:l) 
     integer index2(0:m,0:n,0:l)
     _REAL_   w2(nl,3),w1(ne,nl),w3(nl) &
            ,sd(ne+1),ew(nl),uw(ne,ne),v(nl,nl),w4(ne),ss(nl) 
     integer  ix(nl),iy(nl),iz(nl)
     integer m,n,l,i,j,k
     _REAL_ h,xa,xb,ya,yb,za,zb,coca1,coca2,coca3,cocb1,cocb2,cocb3
     _REAL_ cocc1,cocc2,cocc3,curv,x1,y1,z1,unout,alf
     integer if0,inf,j0,j1,job,k0,k1,nsub,i0,i1,i2
     _REAL_ pi,t,ts
 
      pi = 3.14159265358979d0
      alf=5.1d0*h
      if0 = int(alf/h) + 1

      i0 = nint((x1-xa)/h)
      j0 = nint((y1-ya)/h)
      k0 = nint((z1-za)/h)

      w3 = 0.0d0
      sd = 0.0d0
      ew = 0.0d0
      uw = 0.0d0
      v = 0.0d0
      w4 = 0.0d0

      !------------------------- Generate the matrix --------------------------

      nsub = 0
      do 70 i2=0, if0       !# Starting from the center
         do 10 i1 = i0-i2, i0+i2
            do 10 j1 = j0-i2, j0+i2
               do 10 k1 = k0-i2, k0+i2
                  if(i1<1 .or. j1<1 .or. k1<1) cycle
                  if(i1>m-1 .or. j1>n-1 .or. k1>l-1) cycle
                  if( abs(i1-i0) + abs(j1-j0) + abs(k1-k0) == i2) then
                     t = (x(i1)-x1)*(x(i1)-x1)+(y(j1)-y1)*(y(j1)-y1) &
                           + (z(k1)-z1)*(z(k1)-z1)
                     ts = dsqrt(t)
                     if(ts < alf .and. nsub < nl &
                           .and. phi(i1,j1,k1) > 0.0) then
                        ts = 1.d0+dcos(pi*ts/alf)
                        nsub = nsub + 1
                        ix(nsub) = i1
                        iy(nsub) = j1
                        iz(nsub) = k1
                        ss(nsub) = ts
                        call trans(x1,y1,z1,x(i1),y(j1),z(k1),w2(nsub,1) &
                              ,w2(nsub,2),w2(nsub,3),coca1,coca2,coca3 &
                              ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3)
                        w1(1,nsub) = ts
                        w1(2,nsub) = w2(nsub,1)*ts
                        w1(3,nsub) = w2(nsub,2)*ts
                        w1(4,nsub) = w2(nsub,3)*ts
                        w1(5,nsub) = 0.5d0*w2(nsub,1)*w2(nsub,1)*ts
                        w1(6,nsub) = 0.5d0*w2(nsub,2)*w2(nsub,2)*ts
                        w1(7,nsub) = 0.5d0*w2(nsub,3)*w2(nsub,3)*ts
                        w1(8,nsub) = w2(nsub,1)*w2(nsub,2)*ts
                        w1(9,nsub) = w2(nsub,1)*w2(nsub,3)*ts
                        w1(10,nsub) = w2(nsub,3)*w2(nsub,2)*ts
                     end if
                  end if
         10 continue
      70 continue

      job = 11
      call dsvdc(w1,ne,ne,nsub,sd,ew,uw,ne,v,nl,w4,job,inf)

      do i1=1,ne
         if(abs(sd(i1)) > 1.0d-13) then
            ew(i1)=uw(2,i1)/sd(i1)
         else
            ew(i1)= 0.d0
         end if
      end do

      do i1=1,nl
         w3(i1) = 0.d0
         do j1=1,ne
            w3(i1) = w3(i1) + v(i1,j1)*ew(j1)
         end do
      end do

      unout = 0.d0

      do 80 i2=1, nsub
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)
         unout = unout + w3(i2)*u(i1,j1,k1)*ts
      80 continue

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine interp4i here]
   subroutine interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
         ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout,kdm,kzone)
      implicit none
      integer nl,ne
      parameter( nl = 16, ne = 10 )

      !--------------------------------------------------------------------------
      !-- Subroutine interp3 interpolates u_n^- and u_n^+ of the solution
      !-- of Poisson equation with jumps along interface at (x1,y1) which is a
      !-- point on the interface. The method is the weighted least squares
      !-- interpolation.
      !--
      !-- m,n:          The number of grid in the x- and y-direction.
      !-- x1,y1:        A point on the interface where the interpolation
      !--               take place.
      !-- u(0:m,0:n):   The grid function to be interpolated.
      !--------------------------------------------------------------------------

      

      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l),u(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l),kzone(0:m,0:n,0:l)
      _REAL_ w2(nl,3),w1(ne,nl),w3(nl) &
            ,sd(ne+1),ew(nl),uw(ne,ne),v(nl,nl),w4(ne),ss(nl) 
      integer ix(nl),iy(nl),iz(nl)
      integer m,n,l,i,j,k,kdm
      _REAL_  h,xa,xb,ya,yb,za,zb,coca1,coca2,coca3,cocb1,cocb2,cocb3
      _REAL_ cocc1,cocc2,cocc3,curv,x1,y1,z1,unout,alf
      integer i0,i1,i2,if0,inf,j0,j1,job,k0,k1,nsub
      _REAL_  pi,t,ts    
 
      pi = 3.14159265358979d0
      alf=5.1d0*h
      if0 = int(alf/h) + 1

      i0 = nint((x1-xa)/h)
      j0 = nint((y1-ya)/h)
      k0 = nint((z1-za)/h)

      w3 = 0.0d0
      sd = 0.0d0
      ew = 0.0d0
      uw = 0.0d0
      v = 0.0d0
      w4 = 0.0d0

      !------------------------- Generate the matrix --------------------------

      nsub = 0
      do 70 i2=0, if0       !# Starting from the center
         do 10 i1 = i0-i2, i0+i2
            do 10 j1 = j0-i2, j0+i2
               do 10 k1 = k0-i2, k0+i2
                   if(i1<1 .or. j1<1 .or. k1<1) cycle  
                  if(i1>m-1 .or. j1>n-1 .or. k1>l-1) cycle
                  if( abs(i1-i0) + abs(j1-j0) + abs(k1-k0) == i2) then
                     t = (x(i1)-x1)*(x(i1)-x1)+(y(j1)-y1)*(y(j1)-y1) &
                           + (z(k1)-z1)*(z(k1)-z1)
                     ts = dsqrt(t)
                     if(ts < alf .and. nsub < nl &
                    .and. phi(i1,j1,k1) <=0.0.and.kzone(i1,j1,k1)/100000&
                    .eq. kdm) then
                        ts = 1.d0+dcos(pi*ts/alf)
                        nsub = nsub + 1
                        ix(nsub) = i1
                        iy(nsub) = j1
                        iz(nsub) = k1
                        ss(nsub) = ts
                        call trans(x1,y1,z1,x(i1),y(j1),z(k1),w2(nsub,1) &
                              ,w2(nsub,2),w2(nsub,3),coca1,coca2,coca3 &
                              ,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3)
                        w1(1,nsub) = ts
                        w1(2,nsub) = w2(nsub,1)*ts
                        w1(3,nsub) = w2(nsub,2)*ts
                        w1(4,nsub) = w2(nsub,3)*ts
                        w1(5,nsub) = 0.5d0*w2(nsub,1)*w2(nsub,1)*ts
                        w1(6,nsub) = 0.5d0*w2(nsub,2)*w2(nsub,2)*ts
                        w1(7,nsub) = 0.5d0*w2(nsub,3)*w2(nsub,3)*ts
                        w1(8,nsub) = w2(nsub,1)*w2(nsub,2)*ts
                        w1(9,nsub) = w2(nsub,1)*w2(nsub,3)*ts
                        w1(10,nsub) = w2(nsub,3)*w2(nsub,2)*ts
                     end if
                  end if
         10 continue
      70 continue

      job = 11
      call dsvdc(w1,ne,ne,nsub,sd,ew,uw,ne,v,nl,w4,job,inf)

      do i1=1,ne
         if(abs(sd(i1)) > 1.0d-13) then
            ew(i1)=uw(2,i1)/sd(i1)
         else
            ew(i1)= 0.d0
         end if
      end do

      do i1=1,nl
         w3(i1) = 0.d0
         do j1=1,ne
            w3(i1) = w3(i1) + v(i1,j1)*ew(j1)
         end do
      end do

      unout = 0.d0

      do 80 i2=1, nsub
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)
         unout = unout + w3(i2)*u(i1,j1,k1)*ts
      80 continue

      return
   end subroutine



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine inter_xyi here]
   subroutine inter_xyi(m,n,l,nl,ne,hx,hy,hz,ic,jc,kc &
         ,phi,xp,yp,zp,x,y,z,p,p0,kdm,kzone)

      !# --- This subroutine uses one-sided least squares interpolation to
      !# --- compute p0 of a particular side from a grid p(i,j)
      !# --- Input: m,n,hx,hy,phi,xp,yp,x,y,p,a,b,c,d
      !# --- Output: p0 = p+,

      implicit none
      integer nl,ne
      _REAL_  phi(0:m,0:n,0:l),w1(ne,nl) &
            ,sd(ne+1),ew(nl),uw(ne,ne),vl(nl,nl),w4(ne) &
            ,x(0:m),y(0:n),z(0:l) &
            ,p(0:m,0:n,0:l),w3ux(nl),w3uy(nl),w3uz(nl) &
            ,ewux(nl),ewuy(nl),ewuz(nl) 
       integer ixx(nl),iyy(nl),izz(nl),kzone(0:m,0:n,0:l)
       integer m,n,l,kdm
       _REAL_ hx,hy,hz
       integer ic,jc,kc,i,i1,i2,if0,ii,info,ix,iy,iz,j,job,k,k2,nsub,j2
       _REAL_ xp,yp,zp,p0,alf,h
      h = min(hx,hy,hz)
      alf = 5.1*h
      if0 = int(alf/h)+1
      sd = 0.d0
      ew = 0.d0
      w4 = 0.d0
      w3ux = 0.d0
      w3uy = 0.d0
      w3uz = 0.d0
      ewux = 0.d0
      ewuy = 0.d0
      ewuz = 0.d0
      w1 = 0.d0
      uw = 0.d0
      vl = 0.d0

      !# -------- Interpolates u_{ij} to get u_x^{-} ---------------

     !          ix = nint( (xp - xa)/hx)
     !           iy = nint( (yp - ya)/hy)
     !           iz = nint((zp-za)/hz)
      ix = ic
      iy = jc
      iz = kc

      !# ------------------- Initialize ------------------------------
      nsub = 0

      !# --------------------------------------------------------------

      do  i1=0, if0
         do i2=-i1,i1
            do j2=-i1,i1
               do k2=-i1,i1
                  i = i2 + ix
                  j = j2 + iy
                  k = k2 + iz
      !           if (k<1) then
      !              print *,'xyi error',ix,iy,iz,i2,j2,k2,if0,i1
      !              stop
      !           end if
                  if(i<1.or.j<1.or.k<1) cycle
                  if(i>m-1 .or. j>n-1 .or. k>l-1) cycle
                  if (abs(i2)+abs(j2)+abs(k2) == i1) then

                     !---------- Do coordinates transformation if the grid point is involved ----

                     if(nsub <= nl-1 .and. phi(i,j,k) > 0) then
                        nsub = nsub + 1
                        ixx(nsub) = i
                        iyy(nsub) = j
                        izz(nsub) = k
                        w1(1,nsub) = 1.0d0
                        w1(2,nsub) = x(i) - xp
                        w1(3,nsub) = y(j) - yp
                        w1(4,nsub) = z(k) - zp
                        if(ne .eq. 10) then
                        w1(5,nsub) = 0.5d0*(x(i) - xp)*(x(i) - xp)
                        w1(6,nsub) = 0.5d0*(y(j) - yp)*(y(j) - yp)
                        w1(7,nsub) = 0.5d0*(z(k) - zp)*(z(k) - zp)
                        w1(8,nsub) = (x(i) - xp)*(y(j) - yp)
                        w1(9,nsub) = (x(i) - xp)*(z(k) - zp)
                        w1(10,nsub) = (z(k) - zp)*(y(j) - yp)
                        
                        endif
                     end if
                     !---------- End coordinates transformation -------------------------------
                  end if
               end do
            end do  !  j2=-i1,i1
         end do  !  i2=-i1,i1
      enddo

      do ii =1,nsub
         w3ux(ii) = 0.0
         ewux(ii) = 0.0
      end do
      if(nsub<=24 .and. ne .eq. 10) then
      if(master) then
      write(6,*) 'not enough point for pressure interpolation',nsub,i,j,k
      endif
      endif
      
       
      !----------------- Call least square routine --------------------------

      job = 11
      call dsvdc(w1,ne,ne,nsub,sd,ew,uw,ne,vl,nl,w4,job,info)

      do i=1,ne
         if(abs(sd(i)) > 1.0d-12) then
            ewux(i) = uw(1,i)/sd(i)
         else
            ewux(i) = 0.0d0
         end if
      end do

      do i=1,nsub
         w3ux(i) = 0.0d0
         do j=1,ne
            w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
         end do
      end do

      !# -------- Begin to compute the correction term for u_n^{-} ------
      p0 = 0.d0

      do  i1=1, nsub
         i = ixx(i1)
         j = iyy(i1)
         k = izz(i1)
         p0 = p0 +  w3ux(i1)*p(i,j,k)
      enddo

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a onesubroutine matvec here]
   subroutine matvec(m,n,l,np,n2,n3,istart,indexx,indexy &
         ,indexz,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p &
         ,fw,gw,hw,u1,v1,w1,p1,phi &
         ,cinfox,cinfoy,cinfoz &
         ,index2,cinfo &
         ,fnt,fvec,bf,kzone,nzone,ubc2,aret,fj,gj,hj)

      implicit none 
      _REAL_ aret(100)

  !    ! common /radius/r0
  !    ! common /aret/aret
      !*************************************************************************
      !   This routine
      ! -----------------------------------------------------------------------
      ! x,y,z:coordinate      
      ! u,v,w: velocity
      ! p1,p0,p: pressure
      ! phi: level set
      ! index2: projection point
      ! cinfo: coordinate, normal,tangential of irregular points
      ! ww:working array
      !cinfox,cinfoy,cinfoz:coordinate,normal,tangential of irregular points
      !indexx,indexy,indexz:irregular points
      !ft,fn,fm:jump conditions um,un,ut
      !fec,ubc2:matrix vector multiplication result
      !bf:right-hand-side
      !fnt:merge ft,fn,fm
      !aret:volume of each domain
      !ww1,ww2,vl,uw: SVD
      !kzone: domain of grid points 





      _REAL_ x(0:m),y(0:n),z(0:l)&
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l)&
            ,fw(0:m,0:n,0:l),gw(0:m,0:n,0:l),hw(0:m,0:n,0:l) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l)&
            ,p1(0:m,0:n,0:l),p(0:m,0:n,0:l),p2(0:m,0:n,0:l)&
            ,phi(0:m,0:n,0:l)&
            ,fj(0:m,0:n,0:l),gj(0:m,0:n,0:l),hj(0:m,0:n,0:l) 
     integer index2(0:m,0:n,0:l) 
     _REAL_  cinfo(np*n2) &
            ,cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2) 
     integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
     _REAL_ ft(n2),fn(n2),fm(n2),fvec(n3),ubc2(n2),bf(n3),fnt(n3)
      integer kzone(0:m,0:n,0:l)
      integer i,j,k,l,m,n,np,n2,n3,istart,i00,j00,k00
      _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,rnu,bci1,bci2,bci3
      integer nzone,nlp,nep
      _REAL_ coca1,coca2,cc,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ ct,ct0,curv,ef,h,fuju,fvjv,fwjw
      integer i0,j0,k0,n22,nn,kdm
      _REAL_ pi,pin,pt1,pt2,pt,t1,tf,usout,utout,uvout,vsout,vtout,vvout
      _REAL_ wsout,wtout,wvout,x1,y1,z1,uout,vout,wout 


      !# fnt(n3), the new force combined together
      !read um,un,ut
      do i=1,n2
         ft(i) = fnt(3*i-2)
         fn(i) = fnt(3*i-1)
         fm(i) = fnt(3*i)

      end do
      !calculate velocity
      call bcg(m,n,l,np,n2,n3,istart,indexx,indexy,indexz &
            ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p &
            ,fw,gw,hw,u1,v1,w1,p1,phi,cinfox,cinfoy &
            ,cinfoz,index2,cinfo,ft,fn,fm,fj,gj,hj,kzone)
      !#------ using the coupled interpolating method to evaluate the
      !#------ BC, which leads also the residue for the Schur complement.

      !       call phiout(m,n,gw)
      !       stop
      ! ----------------------------extract coordinate and normal tangential information
      ef = 0.0
      pi = 4.*datan(1.0d0)
      ct = 4.d0*pi/3.d0
      ct0=ct
      if(solverorder .eq. 1) then
      nlp=10
      nep=4
      else
      nlp=27
      nep=10
      endif  
      if(numbubble .eq. 2) then 
      if ( nzone .eq. 1 ) ct = ct*2.0d0
      endif
   !   cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !       print *,"ppp=========="

  !!   do k=0,l
  !!      do j=0,n
  !!         do i=0,m
  !!            p2(i,j,k)=2.0*p(i,j,k)-p1(i,j,k)
  !!         end do
  !!      end do
  !!   end do
   !
      do k=1,l-1
         do j=1,n-1
             do i=1,m-1
                nn = index2(i,j,k)
               n22 = abs(nn)
               ! was only outside, now inside --CQ
               if(nn < 0 ) then
                  !              print *,'matvec',n22,i,j,k,index2(i,j,k)
                  !              x1 = cinfo(n22*(np-1) + 1)
                  !              y1 = cinfo(n22*(np-1) + 2)
                  !              coc = cinfo(n22*(np-1) + 3)
                  !              sic = cinfo(n22*(np-1) + 4)
                  !              curv = cinfo(n22*(np-1) + 5)
                  x1 = cinfo((n22-1)*np  + 1)
                  y1 = cinfo((n22-1)*np  + 2)
                  z1 = cinfo((n22-1)*np  + 3)
                  coca1 = cinfo((n22-1)*np + 4)
                  coca2 = cinfo((n22-1)*np + 5)
                  coca3 = cinfo((n22-1)*np + 6)
                  cocb1 = cinfo((n22-1)*np + 7)
                  cocb2 = cinfo((n22-1)*np + 8)
                  cocb3 = cinfo((n22-1)*np + 9)
                  cocc1 = cinfo((n22-1)*np +10)
                  cocc2 = cinfo((n22-1)*np +11)
                  cocc3 = cinfo((n22-1)*np +12)
                  curv = cinfo( (n22-1)*np +16)
                  kdm = nint(cinfo( (n22-1)*np +17))
                   i00= nint((x1-xa)/hx)
                   j00= nint((y1-ya)/hy)
                   k00= nint((z1-za)/hz)                
                  if(numbubble .eq. 3) then
                  ct=ct0*ctfactor(kdm)
                  endif             
                  !               write(58,'(3i6,3f12.6)') i,j,k,x(i),y(j),z(k)
                  !               write(58,'(3f12.6)') x1,y1,z1
                  !               write(58,'(3f12.6)') coc1,coc2,coc3
                  !--------------------------------pressure on the interface
           call nsetimer_start(4)
                  if(ifupdate .eqv. .true.) then
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i00,j00,k00 &
                        ,phi,x1,y1,z1,x,y,z,p,pt1,kdm,kzone)
                  call inter_xyi(m,n,l,nlp,nep,hx,hy,hz,i00,j00,k00 &
                        ,phi,x1,y1,z1,x,y,z,p1,pt2,kdm,kzone)
                  pt=(1.0d0+ratio_t)*pt1-ratio_t*pt2
                  cinfo(n22*np-np+18)=pt
                  else
                  pt=cinfo(n22*np-np+18)
                  endif
               

           call nsetimer_stop(4)
           call nsetimer_start(5)
                  !          print *,'interp3',n22
                  !------------ut,vt,wt,um,vm,wm,un,vn,wn on the interface
                if(ifupdate .eqv. .true.) then
            call fjmps(fuju,phi,fj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
            call fjmps(fvjv,phi,gj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
            call fjmps(fwjw,phi,hj(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                 cinfo(n22*np-np+19)=fuju
                 cinfo(n22*np-np+20)=fvjv
                 cinfo(n22*np-np+21)=fwjw
                 else   
                 fuju=cinfo(n22*np-np+19)
                 fvjv=cinfo(n22*np-np+20)
                 fwjw=cinfo(n22*np-np+21)
                 endif
                  call interp3d(m,n,l,np,n2,i,j,k,index2,hx,hy,hz &
                        ,xa,xb,ya,yb,za,zb,t,rnu,phi,fw,gw,hw,ft,fn,fm &
                        ,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
                        ,cocc1,cocc2,cocc3,curv,x1,y1,z1 &
                        ,uvout,vvout,wvout,utout,vtout,wtout &
                        ,usout,vsout,wsout,uout,vout,wout,fuju,fvjv,fwjw,kzone)
           call nsetimer_stop(5)
           call nsetimer_start(6)
                  !          call interp3b(m,n,l,np,n2,i,j,k,index2,hx,hy,hz
                  !    1       ,xa,xb,ya,yb,za,zb,t,rnu,phi,fw,gw,hw,ft2,fn2,fm2
                  !    2       ,x,y,z,cinfo,coc1,coc2,coc3,curv,x1,y1,z1
                  !    3       ,utout,vtout,wtout)

                  !          call interp3c(m,n,l,np,n2,i,j,k,index2,hx,hy,hz
                  !    1       ,xa,xb,ya,yb,za,zb,t,rnu,phi,fw,gw,hw,ft2,fn2,fm2
                  !    2       ,x,y,z,cinfo,coc1,coc2,coc3,curv,x1,y1,z1
                  !    3       ,usout,vsout,wsout)

                  !#-------------------right side term for velocity------------------------------------------------

                  !       call gettan(coc1,coc2,coc3,coca1,coca2,coca3
                  !    1    ,cocb1,cocb2,cocb3)
                  !viscosity force on the three directions
                  bci1 = rnu*( uvout*coca1 + vvout*coca2 + wvout*coca3 )
                  bci2 = rnu*( 1.d0*(utout*coca1+vtout*coca2+wtout*coca3) &
                        +1.d0*(uvout*cocb1+vvout*cocb2+wvout*cocb3))
                  bci3 = rnu*( 1.d0*(usout*coca1+vsout*coca2+wsout*coca3) &
                        +1.d0*(uvout*cocc1+vvout*cocc2+wvout*cocc3))

                  !       print *,'coc',coc1,coc2,coc3
                  !       print *,'coc',coca1,coca2,coca3
                  !       print *,'coc',cocb1,cocb2,cocb3
                  !       print *, uvout, vvout, wvout
                  !       print *, 'real',x1,y1,z1
                  !       print *, 'real',-cc*x1/r0**3,-cc*y1/r0**3,-cc*z1/r0**3
                  !       print *, utout, vtout, wtout
                  !       print *, 'real',cc*x1/r0**3,cc*y1/r0**3,cc*z1/r0**3
                  !       print *, usout, vsout, wsout
                  !       print *, 'real',-cc*y1/r0**2/sqrt(x1**2+y1**2)
                  !    1    ,cc*x1/r0**2/sqrt(x1**2+y1**2),0.d0

                  t1 = t+ dt
                  !              fvec(2*n22-1) = uvout - unbc
                  !              fvec(2*n22) = vvout - vnbc

                  !              fvec(2*n22-1) = bci1 - bc1 - 0.0*curv
                  !    1                   + pe(rnu,t,x1,y1)-pt
                  !              fvec(2*n22) = bci2 - bc2

                  pin = ct/(aret(kdm))
                  !              pt = pt + pin
                  !              print *,"ppp",x1,y1,pt,pt+suj,pin
                  !              pt = pt + suj
                  !              pt = pto
                  !              print *,"area",idx,aret(kdm)
                  !              fvec(2*n22-1) = pt-1.0*curv-2.0*bci1
                  !              fvec(2*n22-1) = pt-1.0*curv+2.0*bci1
                  !force balance on the right side
                  fvec(3*n22-2) = pt-1.0d0*curv-2.0d0*bci1

                  ubc2(n22)=pin+1.0d0*curv+2.0d0*bci1
          !        write(444444,*) ubc2(n22),n22
                  if (kzone(i,j,k)<=0) then
                     print *, "errorrrr"
                     stop
                  end if
                  fvec(3*n22-2) = fvec(3*n22-2) - pin
                  fvec(3*n22-1) = bci2
                  fvec(3*n22) = bci3

                  if(nse_debug .eq. 5) then
                  fvec(3*n22-2)=uout-ue(t,rnu,x1,y1,z1)
                  fvec(3*n22-1)=vout-ve(t,rnu,x1,y1,z1)
                  fvec(3*n22)=wout-we(t,rnu,x1,y1,z1)
                  endif
                  
                  if(nse_debug .eq. 6) then
                
                  fvec(3*n22-2)=uvout+2*cc*x1/(r0**4)
                  fvec(3*n22-1)=vvout+2*cc*y1/(r0**4)
                  fvec(3*n22)=wvout+2*cc*z1/(r0**4)
                  endif

                  !              fvec(2*n22-1) = uvout*coc + vvout*sic- unbc*coc
                  !              fvec(2*n22) = vvout - vnbc

                  tf = max( abs(fvec(3*n22-2)), abs(fvec(3*n22-1)) &
                        , abs(fvec(3*n22)))
                  if(tf > ef) then
                     ef = tf
                     i0 = i
                     j0 = j
                     k0=k
                  end if

                  !              print *,'fvec',fvec(3*n22-2),bci2,bci3
                  !              print *,'fvec',curv,pin
                  !              print *,'fvec',-0.5d0*cc*cc/r0**4,pt
                  !              print *,'fvec',-4.d0*rnu*cc/r0**3,2.d0*bci1
                  ! the second term on the right side
                  fvec(3*n22-2) = fvec(3*n22-2) + bf(3*n22-2)
                  fvec(3*n22-1) = fvec(3*n22-1) + bf(3*n22-1)
                  fvec(3*n22) = fvec(3*n22) + bf(3*n22)
               call nsetimer_stop(6)
                  !              write(88,'(f12.6)') fvec(3*n22-2)
                  !              write(88,'(f12.6)') fvec(3*n22-1)
                  !              write(88,'(f12.6)') fvec(3*n22-0)

               end if  ! (nn < 0)
            end do  !  i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1

      return
   end subroutine
   ! make no use

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine gettan here]
   subroutine gettan(coc1,coc2,coc3,coca1,coca2,coca3 &
         ,cocb1,cocb2,cocb3)

      implicit none
      _REAL_ coc1,coc2,coc3,coca1,coca2,coca3,cocb1,cocb2,cocb3
      _REAL_ coca,cocb,cocc,fcta,rmaga      

      if ( abs(abs(coc1) - 1.d0) < 1.d-3 ) then
         coca = 0.d0
         cocb = 1.d0
         cocc = 0.d0
      else
         coca = 1.d0
         cocb = 0.d0
         cocc = 0.d0
      end if

      !       coca = 1.d0
      !       cocb = 0.d0
      !       cocc = 0.d0

      coca1 = coc2*cocc-coc3*cocb
      coca2 = coc3*coca-coc1*cocc
      coca3 = coc1*cocb-coc2*coca
      rmaga = coca1**2+coca2**2+coca3**2
      fcta = dsqrt(1.d0/rmaga)
      coca1 = coca1*fcta
      coca2 = coca2*fcta
      coca3 = coca3*fcta

      cocb1 = coc2*coca3-coc3*coca2
      cocb2 = coc3*coca1-coc1*coca3
      cocb3 = coc1*coca2-coc2*coca1
      !       rmagb = cocb1**2+cocb2**2+cocb3**2
      !       fctb = dsqrt(1.d0/rmagb)
      !       cocb1 = cocb1*fctb
      !       cocb2 = cocb2*fctb
      !       cocb3 = cocb3*fctb

      !       print *, coc1,coc2,coc3
      !       print *, coca1,coca2,coca3
      !       print *, cocb1,cocb2,cocb3
      !       stop

   end subroutine gettan


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine normal here]
   subroutine normal(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1 &
         ,x,y,z,phi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv)
      implicit none

      !-- subroutine normal finds the tangential, normal, curvature information
      !   at a given point on the interface, it may not be a grid point.
      !   Input:  m,n,a,b,c,d,h,x,y,phi
      !           see the description in the subroutine iim_poisson.f
      !      x1,y1:  A point on the interface.
      !   Output: coc,sic,curv,
      !           the normal direction and curvature at (x1,y1)

      _REAL_ phi(0:m,0:n,0:l),x(0:m),y(0:n),z(0:l)
      integer i,j,k,m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1     
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,hxt,hyt,hzt,pha0,phaa0,phab0,phac0
      _REAL_ phb0,phba0,phbb0,phc0,phca0,phcb0,phcc0,phn,phn1,phnc
      _REAL_ phnxy,phnxy1,phnxz,phnxx1,phnxz1,phx0,phx1,phx2
      _REAL_ phbc0,phx3,phx4,phx5,phx6,phx7,phx8,phxx0,phxx1
      _REAL_ phxx2,phxx3,phxx4,phxx5,phxx6,phxx7,phxx8
      _REAL_ phxy0,phxy1,phxy2,phxy3,phxy4,phxy5,phxy6,phxy7,phxy8
      _REAL_ phxz0,phxz1,phxz2,phxz3,phxz4,phxz5,phxz6,phxz7,phxz8
      _REAL_ phy0,phy1,phy2,phy3,phy4,phy5,phy6,phy7,phy8
      _REAL_ phyy0,phyy1,phyy2,phyy3,phyy4,phyy5,phyy6,phyy7,phyy8
      _REAL_ phyz0,phyz1,phyz2,phyz3,phyz4,phyz5,phyz6      
      _REAL_ phyz7,phyz8,phz0,phz1,phz2,phz3,phz4,phz5
      _REAL_ phz6,phz7,phz8,phzz0,phzz1,phzz2,phzz3,phzz4,phzz5
      _REAL_ phzz6,phzz7,phzz8,t1,t2,t3,tem,tp11
      _REAL_ tp12,tp13,tp21,tp22,tp23,tp31,tp32,tp33
      _REAL_ x0,xe,y0,ye,z0,ze

      i = int((x1-xa)/hx)
      j = int((y1-ya)/hy)
      k = int((z1-za)/hz)

      call phiinf(m,n,l,i,j,k,hx,hy,hz,phi,phx1,phy1,phz1 &
            ,phxx1,phyy1,phzz1,phxy1,phxz1,phyz1)
      call phiinf(m,n,l,i+1,j,k,hx,hy,hz,phi,phx2,phy2,phz2 &
            ,phxx2,phyy2,phzz2,phxy2,phxz2,phyz2)
      call phiinf(m,n,l,i,j+1,k,hx,hy,hz,phi,phx3,phy3,phz3 &
            ,phxx3,phyy3,phzz3,phxy3,phxz3,phyz3)
      call phiinf(m,n,l,i+1,j+1,k,hx,hy,hz,phi,phx4,phy4,phz4 &
            ,phxx4,phyy4,phzz4,phxy4,phxz4,phyz4)
      call phiinf(m,n,l,i,j,k+1,hx,hy,hz,phi,phx5,phy5,phz5 &
            ,phxx5,phyy5,phzz5,phxy5,phxz5,phyz5)
      call phiinf(m,n,l,i+1,j,k+1,hx,hy,hz,phi,phx6,phy6,phz6 &
            ,phxx6,phyy6,phzz6,phxy6,phxz6,phyz6)
      call phiinf(m,n,l,i,j+1,k+1,hx,hy,hz,phi,phx7,phy7,phz7 &
            ,phxx7,phyy7,phzz7,phxy7,phxz7,phyz7)
      call phiinf(m,n,l,i+1,j+1,k+1,hx,hy,hz,phi,phx8,phy8,phz8 &
            ,phxx8,phyy8,phzz8,phxy8,phxz8,phyz8)

      !----- Begin nilinear interpolation -------------------------

      hxt = 2.0d0/hx
      hyt = 2.0d0/hy
      hzt = 2.0d0/hz
      x0 = 1.0d0 -( (x1-x(i))*hxt -1.0d0) !# (x1-x(i))*hxt in (-1,1)
      y0 = 1.0d0 -( (y1-y(j))*hyt -1.0d0)
      z0 = 1.0d0 -( (z1-y(k))*hzt -1.0d0)
      xe = 1.0d0 + ( (x1-x(i)) *hxt -1.0d0)
      ye = 1.0d0 + ( (y1-y(j)) *hyt -1.0d0)
      ze = 1.0d0 + ( (z1-z(k)) *hzt -1.0d0)

      phx0 = 0.125*(phx1*x0*y0*z0+phx2*xe*y0*z0+phx3*x0*ye*z0 &
            +phx4*xe*ye*z0+phx5*x0*y0*ze+phx6*xe*y0*ze &
            +phx7*x0*ye*ze+phx8*xe*ye*ze)
      phy0 = 0.125*(phy1*x0*y0*z0+phy2*xe*y0*z0+phy3*x0*ye*z0 &
            +phy4*xe*ye*z0+phy5*x0*y0*ze+phy6*xe*y0*ze &
            +phy7*x0*ye*ze+phy8*xe*ye*ze)
      phz0 = 0.125*(phz1*x0*y0*z0+phz2*xe*y0*z0+phz3*x0*ye*z0 &
            +phz4*xe*ye*z0+phz5*x0*y0*ze+phz6*xe*y0*ze &
            +phz7*x0*ye*ze+phz8*xe*ye*ze)
      phxx0 = 0.125*(phxx1*x0*y0*z0+phxx2*xe*y0*z0+phxx3*x0*ye*z0 &
            +phxx4*xe*ye*z0+phxx5*x0*y0*ze+phxx6*xe*y0*ze &
            +phxx7*x0*ye*ze+phxx8*xe*ye*ze)
      phyy0 = 0.125*(phyy1*x0*y0*z0+phyy2*xe*y0*z0+phyy3*x0*ye*z0 &
            +phyy4*xe*ye*z0+phyy5*x0*y0*ze+phyy6*xe*y0*ze &
            +phyy7*x0*ye*ze+phyy8*xe*ye*ze)
      phzz0 = 0.125*(phzz1*x0*y0*z0+phzz2*xe*y0*z0+phzz3*x0*ye*z0 &
            +phzz4*xe*ye*z0+phzz5*x0*y0*ze+phzz6*xe*y0*ze &
            +phzz7*x0*ye*ze+phzz8*xe*ye*ze)
      phxy0 = 0.125*(phxy1*x0*y0*z0+phxy2*xe*y0*z0+phxy3*x0*ye*z0 &
            +phxy4*xe*ye*z0+phxy5*x0*y0*ze+phxy6*xe*y0*ze &
            +phxy7*x0*ye*ze+phxy8*xe*ye*ze)
      phxz0 = 0.125*(phxz1*x0*y0*z0+phxz2*xe*y0*z0+phxz3*x0*ye*z0 &
            +phxz4*xe*ye*z0+phxz5*x0*y0*ze+phxz6*xe*y0*ze &
            +phxz7*x0*ye*ze+phxz8*xe*ye*ze)
      phyz0 = 0.125*(phyz1*x0*y0*z0+phyz2*xe*y0*z0+phyz3*x0*ye*z0 &
            +phyz4*xe*ye*z0+phyz5*x0*y0*ze+phyz6*xe*y0*ze &
            +phyz7*x0*ye*ze+phyz8*xe*ye*ze)

      !       print *, x1,y1,z1
      !       print *, x0,y0,z0
      !       print *, xe,ye,ze
      !       print *,phx1,phx2,phx3,phx4
      !       print *,phx5,phx6,phx7,phx8

      phn = phx0*phx0 + phy0*phy0 + phz0*phz0
      phn1 = dsqrt(phn) + 1.d-6
      coca1 = phx0/phn1
      coca2 = phy0/phn1
      coca3 = phz0/phn1
      if ( abs(phy0) > abs(phz0) ) then
         phnxy = (phx0*phx0+phy0*phy0)
         phnxy1= dsqrt(phnxy)
         cocb1 = phy0/phnxy1
         cocb2 = -phx0/phnxy1
         cocb3 = 0.d0
         cocc1 = phx0*phz0
         cocc2 = phy0*phz0
         cocc3 = -phnxy
         phnc = dsqrt(cocc1*cocc1+cocc2*cocc2+cocc3*cocc3)
         cocc1 = cocc1/phnc
         cocc2 = cocc2/phnc
         cocc3 = cocc3/phnc
      else
         phnxz = (phx0*phx0+phz0*phz0)
         phnxz1= dsqrt(phnxz)
         cocb1 = phz0/phnxz1
         cocb2 = 0.d0
         cocb3 = -phx0/phnxz1
         cocc1 = -phx0*phy0
         cocc2 = phnxz
         cocc3 = -phy0*phz0
         phnc = dsqrt(cocc1*cocc1+cocc2*cocc2+cocc3*cocc3)
         cocc1 = cocc1/phnc
         cocc2 = cocc2/phnc
         cocc3 = cocc3/phnc
      end if
      !       print *,'local',coca1*cocb1+coca2*cocb2+coca3*cocb3
      !       print *,'local',coca1*cocc1+coca2*cocc2+coca3*cocc3
      !       print *,'local',cocb1*cocc1+cocb2*cocc2+cocb3*cocc3
      !       print *,'local',coca1*coca1+coca2*coca2+coca3*coca3
      !       print *,'local',cocb1*cocb1+cocb2*cocb2+cocb3*cocb3
      !       print *,'local',cocc1*cocc1+cocc2*cocc2+cocc3*cocc3
      !       print *,'local',coca1,coca2,coca3
      !       print *,'local',cocb1,cocb2,cocb3
      !       print *,'local',cocc1,cocc2,cocc3

      call m3v3(1,coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2 &
            ,cocc3,phx0,phy0,phz0,pha0,phb0,phc0)
      call m3n3(1,coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2 &
            ,cocc3,phxx0,phxy0,phxz0,phxy0,phyy0,phyz0,phxz0,phyz0 &
            ,phzz0,tp11,tp12,tp13,tp21,tp22,tp23,tp31,tp32,tp33)
      call m3n3(1,tp11,tp12,tp13,tp21,tp22,tp23,tp31,tp32,tp33 &
            ,coca1,cocb1,cocc1,coca2,cocb2,cocc2,coca3,cocb3,cocc3 &
            ,phaa0,phab0,phac0,phba0,phbb0,phbc0,phca0,phcb0,phcc0)
      !       call m3n3(1,coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2
      !    1    ,cocc3
      !    1    ,coca1,cocb1,cocc1,coca2,cocb2,cocc2,coca3,cocb3,cocc3
      !    1    ,phaa0,phab0,phac0,phba0,phbb0,phbc0,phca0,phcb0,phcc0)
      !       print *, 'local',phaa0,phab0,phac0
      !       print *, 'local',phba0,phbb0,phbc0
      !       print *, 'local',phca0,phcb0,phcc0
      !       stop

      pha0 = pha0 + 1.d-6
      sibb = -phbb0/pha0
      sicc = -phcc0/pha0
      sibc = -phbc0/pha0

      !----- Compute the Curvature --------------------

      t1 = phx0*phx0
      t2 = phy0*phy0
      t3 = phz0*phz0

      ! mean curvature
      curv = (phzz0+phyy0)*t1+(phxx0+phzz0)*t2+(phxx0+phyy0)*t3 &
            - 2.d0*phx0*phy0*phxy0 &
            - 2.d0*phx0*phz0*phxz0 &
            - 2.d0*phz0*phy0*phyz0

      tem = (t1+t2+t3)**1.5d0 + 1.d-6

      curv = -curv/tem/2.d0

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine m3v3 here]
   subroutine m3v3(job,a11,a12,a13,a21,a22,a23,a31,a32,a33 &
         ,v1,v2,v3,p1,p2,p3)

      implicit none
      integer job
      _REAL_ a11,a12,a13,a21,a22,a23,a31,a32,a33,v1,v2,v3,p1,p2,p3

      if ( job == 1 ) then
         ! m*v
         p1 = a11*v1+a12*v2+a13*v3
         p2 = a21*v1+a22*v2+a23*v3
         p3 = a31*v1+a32*v2+a33*v3
      else
         ! v'*m
         p1 = a11*v1+a21*v2+a31*v3
         p2 = a12*v1+a22*v2+a32*v3
         p3 = a13*v1+a23*v2+a33*v3
      end if

   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine m3n3 here]
   subroutine m3n3(job,a11,a12,a13,a21,a22,a23,a31,a32,a33 &
         ,b11,b12,b13,b21,b22,b23,b31,b32,b33 &
         ,p11,p12,p13,p21,p22,p23,p31,p32,p33)

      implicit none
      integer job
      _REAL_ a11,a12,a13,a21,a22,a23,a31,a32,a33 
      _REAL_  b11,b12,b13,b21,b22,b23,b31,b32,b33
      _REAL_   p11,p12,p13,p21,p22,p23,p31,p32,p33

      if ( job == 1 ) then
         ! m*n
         call m3v3(job,a11,a12,a13,a21,a22,a23,a31,a32,a33 &
               ,b11,b21,b31,p11,p21,p31)
         call m3v3(job,a11,a12,a13,a21,a22,a23,a31,a32,a33 &
               ,b12,b22,b32,p12,p22,p32)
         call m3v3(job,a11,a12,a13,a21,a22,a23,a31,a32,a33 &
               ,b13,b23,b33,p13,p23,p33)
      else
         ! n*m
         call m3v3(job,b11,b12,b13,b21,b22,b23,b31,b32,b33 &
               ,a11,a21,a31,p11,p21,p31)
         call m3v3(job,b11,b12,b13,b21,b22,b23,b31,b32,b33 &
               ,a12,a22,a32,p12,p22,p32)
         call m3v3(job,b11,b12,b13,b21,b22,b23,b31,b32,b33 &
               ,b13,b23,b33,p13,p23,p33)
      end if

   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine one_step here]
   subroutine one_step(m,n,l,np,n2,n3,mm,ki,indexx,indexy &
         ,indexz,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p &
         ,fw,gw,hw,u1,v1,w1,p1,phi &
         ,cinfox,cinfoy,cinfoz &
         ,index2,cinfo&
         ,fnt,bf,fvec,kzone,nzone,aret)

      !ki: index of time steps
      implicit none
      _REAL_ aret(100)
!      ! common /radius/ r0


      !x,y,z:coordinate
      !u,v,w:velocity
      !p,p1: pressure
      !phi,index: level set
      !index2,ubc2: projection points
      !cinfo: coordinate,normal,tangential of irregular points
      !cinfox,cinfoy,cinfoz: coordinate,normal,tangential of irregular
      !points
      !indexx,indexy,indexz: irregular points of (x,y,z)
      !fnt,bf,fvec: merge um,un,ut
      !gmat,ipvt,r,hj: coefficient matrix of augmented variables
      !kzone: domain of grid points
      _REAL_  x(0:m),y(0:n),z(0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,p(0:m,0:n,0:l),p1(0:m,0:n,0:l) &
            ,fw(0:m,0:n,0:l),gw(0:m,0:n,0:l),hw(0:m,0:n,0:l) &
            ,fj(0:m,0:n,0:l),gj(0:m,0:n,0:l),hj(0:m,0:n,0:l) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) &
            ,phi(0:m,0:n,0:l)
     integer index(0:m,0:n,0:l) &
            ,index2(0:m,0:n,0:l)
      _REAL_  ubc2(n2) &
            ,cinfo(np*n2) &
            ,cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2) 
      integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
      _REAL_ fnt(n3),bf(n3),fvec(n3),res(n3) 
      integer kzone(0:m,0:n,0:l) 
      integer i,j,k,m,n,l,np,n2,n3,mm,ki,ii,jj
      _REAL_  xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,rnu,e,e1,e2,e_mean,e_abs_mean
      _REAL_ e_mean_irregular,e_mean_regular,e_regular,e_total,e_abs_total
      _REAL_ e_total_regular,elmbda,et,et_regular,hx1,hx2,hy1,hy2,et_x,et_y,et_z
      _REAL_ hz1,hz2,e_relative
      integer nzone,i0,i0_regular,i1,i2,iavg,ierr,ierror,ifirst,ilast,imax
      integer infot,ires,isndcnt,istart,j0,j0_regular,j1,j2,job,k0
      integer k0_regular,k1,k2,lbdcnd,mbdcnd,mdimf,nbdcnd,ndimf
      integer num,num_regular,iter,ictxtA,ictxtB,myrow,mycol,nmaxr,nmaxc
      integer nprow,npcol,mblock,nblock,iarow,iacol,kdm
      _REAL_ pertrb,t1,tol,h,error,err,ct,pi,u_norm 
      _REAL_,allocatable :: gmat(:,:),work(:)    
      integer,allocatable :: ipvt(:),desca(:),descb(:)
#ifdef MPI
        integer ircvcnt(numtasks),idispl(numtasks)
#endif
      character(10) fname
      character(4) str


      ! print n2 (number of irregualr) and n3 (3*n2)
     call nsetimer_start(7)
     call nsetimer_start(10)
   if (master) then
     write(6,*)'one step, n2,n3', n2,n3
     flush(6)
    endif
      istart = ki-1   !if p(t=0) is exact, then istart = ki.
      !compute right-side term
      do i=1,n3
         fnt(i) = 0.0d0                   !# start with zero to get -b
         bf(i) = 0.0d0
      end do
      ifupdate=.true.  
      allocate(fjumpinfo(8*npfjump*n2,16),numfjump(8*n2),interpinfo(npinterp*n2,27),numinterp(n2)) 
      !       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !  ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2 
         !  ,ui,vi,wi,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2
       call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz &
          ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p&
          ,fw,gw,hw &
          ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2,cinfo&
          ,fnt,fvec,bf,kzone,nzone,ubc2,aret,fj,gj,hj)
       ifupdate=.false.
        call nsetimer_stop(10)
        call nsetimer_stop(7) 
            if(ifsetpressure .eqv. .true.) then
               pi=3.1415926535
                ct=4.0d0*pi/3.0d0
               do i=1,n2
                kdm = nint(cinfo( (i-1)*np +17))
               ubc2(i)=ct/aret(kdm)+cinfo((i-1)*np+16)
               enddo
              
      call pressure(m,n,l,np,n2,xa,xb,ya,yb,za,zb &
         ,index2,hx,hy,hz,t,x,y,z,rnu,u,v,w,p &
         ,phi,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,cinfo,ubc2(1:n2),kzone)
                  ifsetpressure = .false.
      do k=0,l
         do j=0,n
            do i=0,m
          p1(i,j,k)=p(i,j,k)
            end do
         end do
      end do
      deallocate(interpinfo,fjumpinfo,numinterp,numfjump)
          return
            endif
   ! store right-side term
     call nsetimer_start(7)
     if(nse_debug .ge. 5) then
      call nsetimer_start(10) 
      do i=1,n3
         bf(i) = -fvec(i)
         !          fnt(i) = bf(i)
         !          fnt(i) = 0.0d0
         !          if (master) write(80,*) bf(i)
      end do
      !first step, form the coefficient matrix
      imax = 3*max(m,n,l)
      tol = 1.0d-6
      if(velocity_solver .eq. 0) then
      ! gmres
      do i=1,n2
      fnt(3*i-2)=0.d0
      fnt(3*i-1)=0.d0
      fnt(3*i)=0.d0
      enddo
     
      call gmres(m,n,l,mm,np,n2,n3,istart,indexx,indexy&
      ,indexz,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p&
      ,fw,gw,hw,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2&
      ,cinfo,imax,fnt,bf,tol&
      ,iter,error,kzone,nzone,ubc2,aret,fj,gj,hj)
      endif
      call nsetimer_stop(10) 
      ! explicit A for augmented variables
#ifdef MPI
     if(velocity_solver .eq. 1) then
     call nsetimer_start(8)
     allocate (desca(9),descb(9))
      ! decompose the matrix to each computer node
    iavg = n3/numtasks
    ires = n3 - iavg*numtasks
    if ( ires == 0 ) then
       ifirst = iavg*mytaskid + 1
       ilast  = iavg*mytaskid + iavg
    else
       if ( (iavg+1)*(mytaskid+1) .le. n3  ) then
          ifirst = (iavg+1)*mytaskid + 1
          ilast  = (iavg+1)*mytaskid + iavg + 1
       else
          ifirst = (iavg+1)*mytaskid+1
          ilast  = n3
       end if
    end if

    nprow=1;npcol=numtasks;
    if(ires ==0) then
    mblock=iavg
    else
    mblock=iavg+1
    endif
    nblock=mblock
     
    call mpi_barrier(mpi_comm_world,ierr)
    call sl_init(ictxtA,nprow,npcol)
    call BLACS_GRIDINFO( ictxtA, nprow, npcol, myrow, mycol )
    nmaxr=numroc(n3,mblock,myrow,0,nprow)
    nmaxc=numroc(n3,nblock,myrow,0,npcol)

    allocate(gmat(nmaxr,nmaxc),ipvt(mblock+nmaxr))
    call descinit( desca,n3,n3, mblock,nblock, 0, 0, ictxtA, max(1,nmaxr),infot )
    call descinit( descb,n3,1, nblock, 1, 0, 0, ictxtA, max(1,nmaxr), infot )
!      call mpi_barrier(mpi_comm_world,ierr)
      !       print *, 'mytaskid',mytaskid,ifirst,ilast
      !       print *, 'mytaskid',(ilast-ifirst+1)*n3,(ifirst-1)*n3+1
!      call mpi_allgather((ilast-ifirst+1)*n3,1,mpi_integer &
!            ,ircvcnt,1,mpi_integer,mpi_comm_world,ierr)
!      call mpi_allgather((ifirst-1)*n3,1,mpi_integer &
!            ,idispl,1,mpi_integer,mpi_comm_world,ierr)
          
      !       if ( mytaskid == 0 ) then
      !         do ii = 1, numtasks
      !           print *, 'count', ircvcnt(ii),idispl(ii)
      !         end do
      !       end if
     ! ifirst = 1; ilast = n3
     !compute the coefficient matrix
     if(ifirst .le. ilast) then
    do j=ifirst,ilast
     do k=1,n3
        fnt(k) = 0.0d0
     end do
     fnt(j) = float(n3)

   
             call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
              ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
              ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2,cinfo&
              ,fnt,fvec,bf,kzone,nzone,ubc2,aret,fj,gj,hj)
!!       ! jf = j
!!       ! write(str,'(i4)') jf
!!      ! str = adjustl(str)
!!       ! write(fname,'(a)') "fvec."//str
!!       ! open(unit=7,file=fname)
      do k=1,n3
    !     gmat(k,j) = fvec(k)
          gmat(k,j-ifirst+1)=fvec(k)/float(n3)
         ! write(7,'(f12.6)') gmat(k,j)
      end do
!!       ! close(7)

    end do
    endif
  
      !collect the result on each node to form the whole matrix
!#ifdef MPI
 ! call mpi_barrier(mpi_comm_world,ierr)
  !   isndcnt = (ilast-ifirst+1)*n3
  !   call mpi_allgatherv(gmat(1,ifirst),isndcnt,mpi_double_precision &
  !         ,gmat(1,1),ircvcnt(1),idispl(1),mpi_double_precision &
  !         ,mpi_comm_world,ierr)
  !   call mpi_barrier(mpi_comm_world,ierr)
     !       call MPI_FINALIZE(ierr)
     !       print *,'gather',mytaskid,ifirst,ilast,isndcnt
!#endif

 
      !100    format(1062f12.6) ! 32*32*32, r=0.8
      !100    format(4878f12.6) ! 64*64*64, r=0.8

   !   if ( master ) then
         ! open(5,file='do.m',status='unknown')
         ! do j=1,n3
         !   write(5,100)(gmat(i,j),i=1,n3)
         ! enddo
         ! close(5)
   !   end if

      ! compute augmented variables
!    call dgefa(gmat,n3,n3,ipvt,infot)

!     ! endif
    do i=1,n3
       ! bf(i) = -fvec(i)
      fnt(i) = bf(i)
    end do
     call nsetimer_stop(8)
 !    job = 0
    call nsetimer_start(9)
  !    do j=1,n3
 !         do i=1,n3
 !         call pdelset(gmatA,i,j,desca,gmat(i,j))
  !     enddo
  !       call pdelset(fntB,j,1,descb,fnt(j))
  !   enddo
      call mpi_barrier(mpi_comm_world,ierr)
   
      call pdgesv(n3,1,gmat,1,1,desca,ipvt,fnt,1,1,descb,infot)
      call mpi_barrier(mpi_comm_world,ierr)
 !     CALL PDLAWRITE( 'scale.dat', n3, 1, fntB, 1, 1, DESCB, 0, 0, work )
 !        open(101,file='scale.dat')
 !        read(101,*) work(1),work(2)
 !        do i=1,n3
 !        read(101,*) fnt(i)
 !        enddo
 !        close(101)
         
    
    deallocate(gmat,ipvt,desca,descb)
      call mpi_barrier(mpi_comm_world,ierr)
    call mpi_bcast(fnt,n3,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
! call dgesl(gmat,n3,n3,ipvt,fnt,job)
   call blacs_gridexit(ictxtA)
   !call blacs_exit(1)
    call nsetimer_stop(9)
    call nsetimer_start(10)

             call matvec(m,n,l,np,n2,n3,istart,indexx,indexy,indexz&
              ,xa,xb,ya,yb,za,zb,hx,hy,hz,dt,t,x,y,z,rnu,u,v,w,p,fw,gw,hw&
              ,u1,v1,w1,p1,phi,cinfox,cinfoy,cinfoz,index2,cinfo&
              ,fnt,fvec,bf,kzone,nzone,ubc2,aret,fj,gj,hj)
      call resid(n3,fvec,bf,res)
      err=dsqrt(prod1(n3,res,res))/dsqrt(prod1(n3,bf,bf))
      if(master) then      
      write(6,*) err,dsqrt(prod1(n3,res,res)),dsqrt(prod1(n3,bf,bf))
      endif
      if(err .gt. 1.0d-2) then
      if(master) then
      write(6,*) "There is deviation in velocity direct solver "
      endif
      endif
    endif
    call nsetimer_stop(10)
#endif
   call nsetimer_start(10)  
    !update u,v,w
     
      do k=0,l
         do j=0,n
            do i=0,m
               u1(i,j,k)=u(i,j,k)
               v1(i,j,k)=v(i,j,k)
               w1(i,j,k)=w(i,j,k)
               u(i,j,k) = fw(i,j,k)
               v(i,j,k) = gw(i,j,k)
               w(i,j,k) = hw(i,j,k)
            end do
         end do
      end do
      !      if(master) print *, "u(24,26,24)",u(24,26,24),u(26,24,24)
      endif
      !Projection step
      !Imcompressible, but now we don't use it0
      !#     update velocity
   if(master) then
   write(6,*) 'finished velocity computation'
   flush(6)
   endif

    
      
      if(nse_debug .le. 8 .and. nse_debug .ne. 0) then
      do k=0,l
         do j=0,n
            do  i=0,m
               index(i,j,k) = 0     !#  regular grid point
               if(phi(i,j,k) == 0.0) index(i,j,k) = 3
               if(phi(i,j,k) > 0.0) then
                  index(i,j,k)=5
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 4
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 4
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 4

               end if

               if(phi(i,j,k) < 0.0) then
                  index(i,j,k)=1
                  if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 2
                  if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 2
                  if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
                        phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 2

               end if
            end do
         end do
      end do
     
      if(nse_debug .ge.5 ) then 
      
      num_regular=0.0
      e_regular=0.0
      e_total_regular=0.0
      num=0.0
      e=0.0
      e_relative=0.0
      e_total=0.0
      et_x=0.0
      et_y=0.0
      et_z=0.0
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               if(index(i,j,k) .ge. 3) then
                  et = sqrt((ue(t,rnu,x(i),y(j),z(k))-u(i,j,k))**2+&
                  (ve(t,rnu,x(i),y(j),z(k))-v(i,j,k))**2+&
                  (we(t,rnu,x(i),y(j),z(k))-w(i,j,k))**2)
                  e_total = abs(et)+e_total

                  num=num+1
                  !           write(53,'(a,3i6,3f12.8)') 'ppp',i,j,k,p(i,j,k),
                  !    1         pe(rnu,t,x(i),y(j),z(k)),et
                  if( abs(et) > abs(e) ) then
                     e = et
                     i0 = i
                     j0 = j
                     k0 = k
                  end if
              if(abs(ue(t,rnu,x(i),y(j),z(k))-u(i,j,k)) .ge. abs(et_x)) then
                et_x= (ue(t,rnu,x(i),y(j),z(k))-u(i,j,k))
               end if

              if(abs(ve(t,rnu,x(i),y(j),z(k))-v(i,j,k)) .ge. abs(et_y)) then
                et_y= (ve(t,rnu,x(i),y(j),z(k))-v(i,j,k))
               end if

              if(abs(we(t,rnu,x(i),y(j),z(k))-w(i,j,k)) .ge. abs(et_z)) then
                et_z= (we(t,rnu,x(i),y(j),z(k))-w(i,j,k))
               end if
               endif
                             
               if(index(i,j,k)==3 .or. index(i,j,k)==4) then
                  et = sqrt((ue(t,rnu,x(i),y(j),z(k))-u(i,j,k))**2+&
                  (ve(t,rnu,x(i),y(j),z(k))-v(i,j,k))**2+&
                  (we(t,rnu,x(i),y(j),z(k))-w(i,j,k))**2)
                  u_norm = sqrt((ue(t,rnu,x(i),y(j),z(k)))**2+&
                  (ve(t,rnu,x(i),y(j),z(k)))**2+&
                  (we(t,rnu,x(i),y(j),z(k)))**2)
               if(e_relative .le. et/u_norm) then
               e_relative=et/u_norm 
               endif 
               endif

               if(index(i,j,k) == 5) then
                  !              vi(i,j,k) =  pe(rnu,t,x(i),y(j),z(k))-p(i,j,k)
                  et_regular = abs(ue(t,rnu,x(i),y(j),z(k))-u(i,j,k))+&
                  abs(ve(t,rnu,x(i),y(j),z(k))-v(i,j,k))+&
                  abs(we(t,rnu,x(i),y(j),z(k))-w(i,j,k))
                  e_total_regular = et_regular &
                        +e_total_regular
                  num_regular=num_regular+1
                  !           write(53,'(a,3i6,3f12.8)') 'ppp',i,j,k,p(i,j,k),
                  !    1         pe(rnu,t,x(i),y(j),z(k)),et
                  if( abs(et_regular) > abs(e_regular)) then
                     e_regular = et_regular
                     i0_regular = i
                     j0_regular = j
                     k0_regular = k
                  end if
                    
               end if
            end do  !  i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1
      e_mean=e_total/num
      e_mean_regular=e_total_regular/num_regular
      e_mean_irregular=(e_total-e_total_regular)/(num-num_regular)
      if(master) then
       write(6,*)&
            "checking velocity give analytical pressure and velocity"
       write(6,"(A,ES14.2E2)")&
            'max error component in velocity=  ',abs(et_x)+abs(et_y)+abs(et_z)&
            ,'max error in velocity= ',e&
            ,'relative error in velocity',e_relative& 
            ,'mean error in velocity=  ',e_mean
       write(6,"(A,ES14.4E2)")&
            'component of x=  ' ,ue(t,rnu,x(i0),y(j0),z(k0)) &
            ,'component of y=  ' ,ve(t,rnu,x(i0),y(j0),z(k0)) &
            ,'component of z=  ',we(t,rnu,x(i0),y(j0),z(k0)) 
      write(6,*)&
            'coordinate=  ',i0,j0,k0
      endif
      if(nse_debug .le. 7) then
      return
      endif
      endif
      endif
   
  !update pressure for k-1 step
      do k=0,l
         do j=0,m
            do i=0,n
               p1(i,j,k)=p(i,j,k)
            end do
         end do
      end do
     call nsetimer_stop(10) 
  !    do i=1,n2
  !       write(666666,*) ubc2(i)
  !    end do
          
      ! ui, vi, wi are initialized in pressure --CQ
      !compute pressure
      call pressure(m,n,l,np,n2,xa,xb,ya,yb,za,zb,index2 &
            ,hx,hy,hz,t,x,y,z,rnu,u,v,w,p,phi &
            ,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,cinfo,ubc2,kzone)
      call mpi_barrier(mpi_comm_world,ierr)
      call nsetimer_start(10)
      deallocate(interpinfo,fjumpinfo,numfjump,numinterp)
      if(nse_debug .le.3 .and. nse_debug .ne. 0) then
          return
       endif
         !check pressure

     if(nse_debug .eq. 4 .or. nse_debug .eq. 8)then 

      num_regular=0.0
      e_regular=0.0
      e_total_regular=0.0
      num=0.0
      e_relative=0.0  
      e=0.0
      e_total=0.0
      e_abs_total=0.0
      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               ! only outside --CQ
               if(index(i,j,k) .ge. 3) then
                  !              vi(i,j,k) =  pe(rnu,t,x(i),y(j),z(k))-p(i,j,k)
                  et = abs(pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))
                  e_total = (pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))+e_total
                  e_abs_total =abs (pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))+e_abs_total
                    
                  num=num+1
                  !           write(53,'(a,3i6,3f12.8)') 'ppp',i,j,k,p(i,j,k),
                  !    1         pe(rnu,t,x(i),y(j),z(k)),et
                  if( abs(et) > abs(e) ) then
                     e = et
                     i0 = i
                     j0 = j
                     k0 = k
                  end if
               end if
               if(index(i,j,k)==3 .or. index(i,j,k)==4) then
               if(e_relative .le. abs((pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))/pe(rnu,t,x(i),y(j),z(k)))) then
               e_relative=abs((pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))/pe(rnu,t,x(i),y(j),z(k)))
               endif 
               endif


               if(index(i,j,k) == 5) then
                  !              vi(i,j,k) =  pe(rnu,t,x(i),y(j),z(k))-p(i,j,k)
                  et_regular = -(pe(rnu,t,x(i),y(j),z(k))-p(i,j,k))
                  e_total_regular = abs(pe(rnu,t,x(i),y(j),z(k))-p(i,j,k)) &
                        +e_total_regular
                  num_regular=num_regular+1
                  !           write(53,'(a,3i6,3f12.8)') 'ppp',i,j,k,p(i,j,k),
                  !    1         pe(rnu,t,x(i),y(j),z(k)),et
                  if( abs(et_regular) > abs(e_regular)) then
                     e_regular = et_regular
                     i0_regular = i
                     j0_regular = j
                     k0_regular = k
                  end if
               end if
            end do  !  i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1
      e_mean=e_total/num
      e_abs_mean=e_abs_total/num
      e_mean_regular=e_total_regular/num_regular
      e_mean_irregular=(e_total-e_total_regular)/(num-num_regular)
       if(master) then
       if(nse_debug .eq.4) then
       write(6,*)&
           'checking pressure after gmres with analytical sources &
                   and interface pressure'  
       else
       write(6,*)&
            'checking pressure'
       endif
       write(6,'(A,ES14.2E2)') 'max error in pressure=  ', e&  
            ,'relative error  ',e_relative
       write(6,'(A,ES14.4E2)')&
       'analytical pressure =',pe(rnu,t,x(i0),y(j0),z(k0)), &
            'real pressure =',p(i0,j0,k0)
        write(6,*)'coordinate= ',i0,j0,k0
        write(6,'(A,ES14.2E2)')&
       'mean  error= ',e_abs_mean,&
       'mean system error = ',e_mean
       flush(6) 
       endif
    endif
      call mpi_barrier(mpi_comm_world,ierr)
     call nsetimer_stop(10)
     call nsetimer_stop(7)
      

      
   end subroutine one_step



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiinf here]
   subroutine phiinf(m,n,l,i,j,k,hx,hy,hz,phi,phx,phy,phz &
         ,phxx,phyy,phzz,phxy,phxz,phyz)
      implicit none

      !#-- subroutine phiinf find the derivatives information of a lvel set
      !#   phi(x,y) at a grid point (x_i,y_j) using the central difference.

      _REAL_ phi(0:m,0:n,0:l)
      integer m,n,l,i,j,k
      _REAL_ hx,hy,hz,phx,phy,phz,phxx,phyy,phzz,phxy,phxz,phyz     
      _REAL_ h1xy,h1xz,h1yz,hx1,hx2,hy1,hy2,hz1,hz2 
      

      hx2 = 2.0d0*hx
      hx1 = hx*hx
      hy2 = 2.0d0*hy
      hy1 = hy*hy
      hz2 = 2.0d0*hz
      hz1 = hz*hz
      h1xy = hx*hy
      h1xz = hx*hz
      h1yz = hz*hy

      if( i == 0 ) then
         phx = (-3.d0*phi(i,j,k)+4.d0*phi(i+1,j,k)-phi(i+2,j,k))/hx2
         phxx=(phi(i+2,j,k)+phi(i,j,k)-2.d0*phi(i+1,j,k))/hx1
      else
         if(i == m) then
            phx = (3.d0*phi(i,j,k)-4.d0*phi(i-1,j,k)+phi(i-2,j,k))/hx2
            phxx=(phi(i-2,j,k)+phi(i,j,k)-2.d0*phi(i-1,j,k))/hx1
         else
            phx = (phi(i+1,j,k) - phi(i-1,j,k))/hx2
            phxx=(phi(i+1,j,k)+phi(i-1,j,k)-2.d0*phi(i,j,k))/hx1
         end if
      end if

      if( j == 0 ) then
         phy = (-3.d0*phi(i,j,k)+4.d0*phi(i,j+1,k)-phi(i,j+2,k))/hy2
         phyy=(phi(i,j+2,k)+phi(i,j,k)-2.d0*phi(i,j+1,k))/hy1
      else
         if(j == n) then
            phy= (3.d0*phi(i,j,k)-4.d0*phi(i,j-1,k)+phi(i,j-2,k))/hy2
            phyy=(phi(i,j-2,k)+phi(i,j,k)-2.d0*phi(i,j-1,k))/hy1
         else
            phy = (phi(i,j+1,k) - phi(i,j-1,k))/hy2
            phyy= (phi(i,j+1,k)+phi(i,j-1,k)-2.d0*phi(i,j,k))/hy1
         end if
      end if

      if( k == 0 ) then
         phz = (-3.d0*phi(i,j,k)+4.d0*phi(i,j,k+1)-phi(i,j,k+2))/hz2
         phzz=(phi(i,j,k+2)+phi(i,j,k)-2.d0*phi(i,j,k+1))/hz1
      else
         if(k == l) then
            phz= (3.d0*phi(i,j,k)-4.d0*phi(i,j,k-1)+phi(i,j,k-2))/hz2
            phzz=(phi(i,j,k-2)+phi(i,j,k)-2.d0*phi(i,j,k-1))/hz1
         else
            phz = (phi(i,j,k+1) - phi(i,j,k-1))/hz2
            phzz= (phi(i,j,k+1)+phi(i,j,k-1)-2.d0*phi(i,j,k))/hz1
         end if
      end if

      if( i == 0 .and. j == 0) then
         phxy = (phi(i+2,j+2,k)+phi(i,j,k)-phi(i+2,j,k) &
               -phi(i,j+2,k))/(4.0d0*h1xy)
      else if( i == 0 .and. j == n) then
         phxy = (phi(i,j-2,k)+phi(i+2,j,k)-phi(i+2,j-2,k) &
               -phi(i,j,k))/(4.0d0*h1xy)
      else if( i == m .and. j == 0) then
         phxy = (phi(i-2,j,k)+phi(i,j+2,k)-phi(i-2,j+2,k) &
               -phi(i,j,k))/(4.0d0*h1xy)
      else if( i == m .and. j == n) then
         phxy = -(phi(i,j-2,k)+phi(i-2,j,k)-phi(i-2,j-2,k) &
               -phi(i,j,k))/(4.0d0*h1xy)
      else if (i == 0) then
         phxy = (phi(i+2,j+1,k)+phi(i,j-1,k)-phi(i+2,j-1,k) &
               -phi(i,j+1,k))/(4.0d0*h1xy)
      else if (i == m) then
         phxy = (phi(i,j+1,k)+phi(i-2,j-1,k)-phi(i,j-1,k) &
               -phi(i-2,j+1,k))/(4.0d0*h1xy)
      else if (j == 0) then
         phxy = (phi(i+1,j+2,k)+phi(i-1,j,k)-phi(i+1,j,k) &
               -phi(i-1,j+2,k))/(4.0d0*h1xy)
      else if (j == n) then
         phxy = (phi(i+1,j,k)+phi(i-1,j-2,k)-phi(i+1,j-2,k) &
               -phi(i-1,j,k))/(4.0d0*h1xy)
      else
         phxy = (phi(i+1,j+1,k)+phi(i-1,j-1,k)-phi(i+1,j-1,k) &
               -phi(i-1,j+1,k))/(4.0d0*h1xy)
      end if

      if( i == 0 .and. k == 0) then
         phxz = (phi(i+2,j,k+2)+phi(i,j,k)-phi(i+2,j,k) &
               -phi(i,j,k+2))/(4.0d0*h1xz)
      else if( i == 0 .and. k == l) then
         phxz = (phi(i,j,k-2)+phi(i+2,j,k)-phi(i+2,j,k-2) &
               -phi(i,j,k))/(4.0d0*h1xz)
      else if( i == m .and. k == 0) then
         phxz = (phi(i-2,j,k)+phi(i,j,k+2)-phi(i-2,j,k+2) &
               -phi(i,j,k))/(4.0d0*h1xz)
      else if( i == m .and. k == l) then
         phxz = -(phi(i,j,k-2)+phi(i-2,j,k)-phi(i-2,j,k-2) &
               -phi(i,j,k))/(4.0d0*h1xz)
      else if (i == 0) then
         phxz = (phi(i+2,j,k+1)+phi(i,j,k-1)-phi(i+2,j,k-1) &
               -phi(i,j,k+1))/(4.0d0*h1xz)
      else if (i == m) then
         phxz = (phi(i,j,k+1)+phi(i-2,j,k-1)-phi(i,j,k-1) &
               -phi(i-2,j,k+1))/(4.0d0*h1xz)
      else if (k == 0) then
         phxz = (phi(i+1,j,k+2)+phi(i-1,j,k)-phi(i+1,j,k) &
               -phi(i-1,j,k+2))/(4.0d0*h1xz)
      else if (k == l) then
         phxz = (phi(i+1,j,k)+phi(i-1,j,k-2)-phi(i+1,j,k-2) &
               -phi(i-1,j,k))/(4.0d0*h1xz)
      else
         phxz = (phi(i+1,j,k+1)+phi(i-1,j,k-1)-phi(i+1,j,k-1) &
               -phi(i-1,j,k+1))/(4.0d0*h1xz)
      end if

      if( k == 0 .and. j == 0) then
         phyz = (phi(i,j+2,k+2)+phi(i,j,k)-phi(i,j,k+2) &
               -phi(i,j+2,k))/(4.0d0*h1yz)
      else if( k == 0 .and. j == n) then
         phyz = (phi(i,j-2,k)+phi(i,j,k+2)-phi(i,j-2,k+2) &
               -phi(i,j,k))/(4.0d0*h1yz)
      else if( k == l .and. j == 0) then
         phyz = (phi(i,j,k-2)+phi(i,j+2,k)-phi(i,j+2,k-2) &
               -phi(i,j,k))/(4.0d0*h1yz)
      else if( k == l .and. j == n) then
         phyz = -(phi(i,j-2,k)+phi(i,j,k-2)-phi(i,j-2,k-2) &
               -phi(i,j,k))/(4.0d0*h1yz)
      else if (k == 0) then
         phyz = (phi(i,j+1,k+2)+phi(i,j-1,k)-phi(i,j-1,k+2) &
               -phi(i,j+1,k))/(4.0d0*h1yz)
      else if (k == l) then
         phyz = (phi(i,j+1,k)+phi(i,j-1,k-2)-phi(i,j-1,k) &
               -phi(i,j+1,k-2))/(4.0d0*h1yz)
      else if (j == 0) then
         phyz = (phi(i,j+2,k+1)+phi(i,j,k-1)-phi(i,j,k+1) &
               -phi(i,j+2,k-1))/(4.0d0*h1yz)
      else if (j == n) then
         phyz = (phi(i,j,k+1)+phi(i,j-2,k-1)-phi(i,j-2,k+1) &
               -phi(i,j,k-1))/(4.0d0*h1yz)
      else
         phyz = (phi(i,j+1,k+1)+phi(i,j-1,k-1)-phi(i,j-1,k+1) &
               -phi(i,j+1,k-1))/(4.0d0*h1yz)
      end if

      return
   end subroutine

   !***********************************************************************
   !                          SUBROUTINE PHIINF1
   !***********************************************************************
   !--

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiinf1 here]
   subroutine phiinf1(m,n,l,i,j,k,hx,hy,hz,phi,phx,phy,phz)
      implicit none

      !-- subroutine phiinf find the first derivatives information of a lvel set
      !   phi(x,y) at a grid point (x_i,y_j) using the central deifference.

     _REAL_ phi(0:m,0:n,0:l)
     integer m,n,l,i,j,k
     _REAL_ hx,hy,hz,phx,phy,phz,hx2,hy2,hz2



      hx2 = 2.0*hx
      hy2 = 2.0*hy
      hz2 = 2.0*hz

      if( i == 0 ) then
         phx = (-3.d0*phi(i,j,k)+4.d0*phi(i+1,j,k)-phi(i+2,j,k))/hx2
      else
         if(i == m) then
            phx = (3.d0*phi(i,j,k)-4.d0*phi(i-1,j,k)+phi(i-2,j,k))/hx2
         else
            ! phx = (phi(i+1,j,k) - phi(i-1,j,k))/hx2
            if (phi(i+1,j,k)+phi(i-1,j,k) > phi(i,j,k)+phi(i,j,k)) then
               phx = (phi(i+1,j,k) - phi(i,j,k))/hx
            else
               phx = (phi(i,j,k) - phi(i-1,j,k))/hx
            end if
         end if
      end if

      if( j == 0 ) then
         phy = (-3.d0*phi(i,j,k)+4.d0*phi(i,j+1,k)-phi(i,j+2,k))/hy2
      else
         if(j == n) then
            phy= (3.d0*phi(i,j,k)-4.d0*phi(i,j-1,k)+phi(i,j-2,k))/hy2
         else
            ! phy = (phi(i,j+1,k) - phi(i,j-1,k))/hy2
            if (phi(i,j+1,k)+phi(i,j-1,k) > phi(i,j,k)+phi(i,j,k)) then
               phy = (phi(i,j+1,k) - phi(i,j,k))/hy
            else
               phy = (phi(i,j,k) - phi(i,j-1,k))/hy
            end if
         end if
      end if

      if( k == 0 ) then
         phz = (-3.d0*phi(i,j,k)+4.d0*phi(i,j,k+1)-phi(i,j,k+2))/hz2
      else
         if(k == l) then
            phz= (3.d0*phi(i,j,k)-4.d0*phi(i,j,k-1)+phi(i,j,k-2))/hz2
         else
            ! phz = (phi(i,j,k+1) - phi(i,j,k-1))/hz2
            if (phi(i,j,k+1)+phi(i,j,k-1) > phi(i,j,k)+phi(i,j,k)) then
               phz = (phi(i,j,k+1) - phi(i,j,k))/hz
            else
               phz = (phi(i,j,k) - phi(i,j,k-1))/hz
            end if
         end if
      end if
      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pressure here]t

   subroutine pressure(m,n,l,np,n2,xa,xb,ya,yb,za,zb &
         ,index2,hx,hy,hz,t,x,y,z,rnu,u,v,w,p &
         ,phi,indexx,indexy,indexz,cinfox,cinfoy,cinfoz,cinfo,ubc2,kzone)

      implicit none
     
      !x,y,z:coordinate/
      !u,v,w: velocity
      !p:pressure
      !index2:projection points
      !phi:level set
      !bda,bdb,bdc,bdd,bde,bdf:boundary condition of velocity
      !indexx,indexy,indexz:irregular points
      !cinfox,cinfoy,cinfoz:informaion on irregular points


      _REAL_  x(0:m),y(0:n),z(0:l),&
            u(0:m,0:n,0:l),v(0:m,0:n,0:l) &
            ,w(0:m,0:n,0:l) ,p(0:m,0:n,0:l),bv(0:m,0:n,0:l) 
       integer index2(0:m,0:n,0:l),index(0:m,0:n,0:l),kzone(0:m,0:n,0:l) 
       _REAL_ phi(0:m,0:n,0:l) 
        integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
        _REAL_ cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2),ubc2(n2),cinfo(np*n2)
        integer m,n,l,np,n2,i,i0,j,j0,k,k0,mm2,nx,ny,nz,kdm      
        _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,t,rnu,elmbda,h,hx1,hx2,cc,x1,y1,z1,r1
        _REAL_ a1,b1,coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv
        _REAL_ unout,unin,vnout,vnin,wnout,wnin,tt,wt,xt,yt,zt
        _REAL_ hy1,hy2,hz1,hz2,rt,tux,tuy,tuz,tvx,tvy,tvz,twx,twy,twz,err
        _REAL_,allocatable :: gmat(:,:),bf(:),fn(:),fvec(:),res(:)   
      integer,allocatable :: ipvt(:),desca(:),descb(:)
        integer ierr,ictxtA,infot,iavg,ires 
        integer ifirst,ilast,mblock,nblock,mycol,myrow,nmaxr,nmaxc,nprow,npcol     

      !setting step size
      call nsetimer_start(10)
      h = min(hx,hy,hz)
      hx1 = hx*hx
      hx2 = 2.0d0*hx
      hy1 = hy*hy
      hy2 = 2.0d0*hy
      hz1 = hz*hz
      hz2 = 2.0d0*hz

      !       r0 = 0.8d0
      !       rnu = 0.5d0
      !       cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)

      !      t1 = t + hx

      do k=0,l
         do j=0,n
            do i=0,m
               bv(i,j,k) = 0.0d0
               !           fw(i,j) = u(i,j)
               !           gw(i,j) = v(i,j)
            end do
         end do
      end do

      do  k= 1,l-1
         do  j= 1,n-1
            do  i= 1,m-1

               !           if( phi(i,j)  .gt. 0.0) then   !# if 1
               !# -------------------- Begin x-direction ----------------
              if(phi(i,j,k) .le. 0) then
              bv(i,j,k)=0.0
              else

               if(phi(i,j,k)*phi(i-1,j,k) <= 0.0) then
                  if(phi(i,j,k)*phi(i+1,j,k)>0.0) then 
                  tux = -(u(i+1,j,k)-u(i,j,k))/hx
                  tvx = -(v(i+1,j,k)-v(i,j,k))/hx
                  twx = -(w(i+1,j,k)-w(i,j,k))/hx
                  else
                  tux=0.0d0
                  tvx=0.0d0
                  twx=0.0d0
                  endif
               
       !             nx = indexx(i-1,j,k)
       !             if (phi(i-2,j,k) > 0.0) nx = nx+1
       !             x1 = cinfox(nx*np-np+1)
       !             y1 = cinfox(nx*np-np+2)
       !             z1 = cinfox(nx*np-np+3)
       !             coca1 = cinfox(nx*np-np+4)
       !             coca2 = cinfox(nx*np-np+5)
       !             coca3 = cinfox(nx*np-np+6)
       !             cocb1 = cinfox(nx*np-np+7)
       !             cocb2 = cinfox(nx*np-np+8)
       !             cocb3 = cinfox(nx*np-np+9)
       !             cocc1 = cinfox(nx*np-np+10)
       !             cocc2 = cinfox(nx*np-np+11)
       !             cocc3 = cinfox(nx*np-np+12)
       !             kdm=nint(cinfox(nx*np-np+22))
       !             r1 = x1*x1 + y1*y1 + z1*z1
       !             xt = x1/(r1*r1*r1)
       !             yt = y1/(r1*r1*r1)
       !             zt = z1/(r1*r1*r1)
       !             
       !             ! calculate the jump condition of u
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
       !             wt = (unout - unin)*abs(coca1)
       !             a1 = (x1-x(i-1))/hx
       !             b1 = 1.d0 - a1
       !             tux =-( (u(i,j,k)-u(i-1,j,k)-tt)/hx+a1*wt)

       !            
       !             ! calculate the jump condition of v
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
       !             wt = (vnout - vnin)*abs(coca1)
       !             a1 = (x1-x(i-1))/hx
       !             b1 = 1.d0 - a1
       !             tvx = -((v(i,j,k)-v(i-1,j,k)-tt)/hx+a1*wt)


       !             
       !            ! calculate the jump condition of w
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
       !             wt = (wnout - wnin)*abs(coca1)
       !             a1 = (x1-x(i-1))/hx
       !             b1 = 1.d0 - a1
       !             twx = -((w(i,j,k)-w(i-1,j,k)-tt)/hx+a1*wt)


       !          endif 
               else if(phi(i,j,k)*phi(i+1,j,k) <= 0.0) then     !# if 2
                  tux = -(u(i,j,k)-u(i-1,j,k))/hx
                  tvx = -(v(i,j,k)-v(i-1,j,k))/hx
                  twx = -(w(i,j,k)-w(i-1,j,k))/hx
                  !              elseif(phi(i,j)*phi(i,j-1) .le. 0.0 .or.
                  !    1                phi(i,j)*phi(i,j+1) .le. 0.0) then     !# if 2
               else
                  tux = -(u(i+1,j,k)-u(i-1,j,k))/hx2
                  tvx = -(v(i+1,j,k)-v(i-1,j,k))/hx2
                  twx = -(w(i+1,j,k)-w(i-1,j,k))/hx2
               end if                                          !# if 2

               !#  done with x-direction, now in y-direction  ------------------

               if(phi(i,j,k)*phi(i,j-1,k) <= 0.0) then
                  if(phi(i,j,k)*phi(i,j+1,k)>0.0) then 
                  tuy = -(u(i,j+1,k)-u(i,j,k))/hy
                  tvy = -(v(i,j+1,k)-v(i,j,k))/hy
                  twy = -(w(i,j+1,k)-w(i,j,k))/hy
                  else
                  tuy=0.0d0
                  tvy=0.0d0
                  twy=0.0d0
                  endif
        
       !             ny = indexy(i,j-1,k)
       !             if (phi(i,j-2,k) > 0.0) ny = ny+1
       !             x1 = cinfoy(ny*np-np+1)
       !             y1 = cinfoy(ny*np-np+2)
       !             z1 = cinfoy(ny*np-np+3)
       !             coca1 = cinfoy(ny*np-np+4)
       !             coca2 = cinfoy(ny*np-np+5)
       !             coca3 = cinfoy(ny*np-np+6)
       !             cocb1 = cinfoy(ny*np-np+7)
       !             cocb2 = cinfoy(ny*np-np+8)
       !             cocb3 = cinfoy(ny*np-np+9)
       !             cocc1 = cinfoy(ny*np-np+10)
       !             cocc2 = cinfoy(ny*np-np+11)
       !             cocc3 = cinfoy(ny*np-np+12)
       !             kdm=nint(cinfoy(ny*np-np+22))
       !             r1 = x1*x1 + y1*y1 + z1*z1
       !             xt = x1/(r1*r1*r1)
       !             yt = y1/(r1*r1*r1)
       !             zt = z1/(r1*r1*r1)
       !             ! calculate the jump condition of u
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
       !             wt = (unout - unin)*abs(coca2)
       !             a1 = (y1-y(j-1))/hy
       !             b1 = 1.d0 - a1
       !             tuy =-( (u(i,j,k)-u(i,j-1,k)-tt)/hy+a1*wt)

       !            
       !             ! calculate the jump condition of v
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
       !             wt = (vnout - vnin)*abs(coca2)
       !             a1 = (y1-y(j-1))/hy
       !             b1 = 1.d0 - a1
       !             tvy = -((v(i,j,k)-v(i,j-1,k)-tt)/hy+a1*wt)


       !             
       !            ! calculate the jump condition of w
       !             call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
       !             call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
       !                   ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
       !                   ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
       !             wt = (wnout - wnin)*abs(coca2)
       !             a1 = (y1-y(j-1))/hy
       !             b1 = 1.d0 - a1
       !             twy = -((w(i,j,k)-w(i,j-1,k)-tt)/hy+a1*wt)
       !          endif
               else if(phi(i,j,k)*phi(i,j+1,k) <= 0.0) then
                  tuy = -(u(i,j,k)-u(i,j-1,k))/hy
                  tvy = -(v(i,j,k)-v(i,j-1,k))/hy
                  twy = -(w(i,j,k)-w(i,j-1,k))/hy
                  !              elseif(phi(i,j)*phi(i-1,j) .le. 0.0 .or.
                  !    1                phi(i,j)*phi(i+1,j) .le. 0.0) then
               else
                  tuy = -(u(i,j+1,k)-u(i,j-1,k))/hy2
                  tvy = -(v(i,j+1,k)-v(i,j-1,k))/hy2
                  twy = -(w(i,j+1,k)-w(i,j-1,k))/hy2
               end if

               !#  done with y-direction, now in z-direction  ------------------

               if(phi(i,j,k)*phi(i,j,k-1) <= 0.0) then
                  if(phi(i,j,k)*phi(i,j,k+1) > 0.0) then
                  tuz = -(u(i,j,k+1)-u(i,j,k))/hz
                  tvz = -(v(i,j,k+1)-v(i,j,k))/hz
                  twz = -(w(i,j,k+1)-w(i,j,k))/hz
                  else
                  tuz=0.0d0
                  tvz=0.0d0
                  twz=0.0d0
                  endif

      !              nz = indexz(i,j,k-1)
      !              if (phi(i,j,k-2) > 0.0) nz = nz+1
      !              x1 = cinfoz(nz*np-np+1)
      !              y1 = cinfoz(nz*np-np+2)
      !              z1 = cinfoz(nz*np-np+3)
      !              coca1 = cinfoz(nz*np-np+4)
      !              coca2 = cinfoz(nz*np-np+5)
      !              coca3 = cinfoz(nz*np-np+6)
      !              cocb1 = cinfoz(nz*np-np+7)
      !              cocb2 = cinfoz(nz*np-np+8)
      !              cocb3 = cinfoz(nz*np-np+9)
      !              cocc1 = cinfoz(nz*np-np+10)
      !              cocc2 = cinfoz(nz*np-np+11)
      !              cocc3 = cinfoz(nz*np-np+12)
      !              kdm=nint(cinfoz(nz*np-np+22))
      !              r1 = x1*x1 + y1*y1 + z1*z1
      !              xt = x1/(r1*r1*r1)
      !              yt = y1/(r1*r1*r1)
      !              zt = z1/(r1*r1*r1)
      !              ! calculate the jump condition of u
      !              call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unout)
      !              call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,u,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,unin,kdm,kzone)
      !              wt = (unout - unin)*abs(coca3)
      !              a1 = (z1-z(k-1))/hz
      !              b1 = 1.d0 - a1
      !              tuz =-( (u(i,j,k)-u(i,j,k-1)-tt)/hz+a1*wt)

      !             
      !              ! calculate the jump condition of v
      !              call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnout)
      !              call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,v,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,vnin,kdm,kzone)
      !              wt = (vnout - vnin)*abs(coca3)
      !              a1 = (z1-z(k-1))/hz
      !              b1 = 1.d0 - a1
      !              tvz = -((v(i,j,k)-v(i,j,k-1)-tt)/hz+a1*wt)


      !              
      !             ! calculate the jump condition of w
      !              call interp4(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnout)
      !              call interp4i(m,n,l,i,j,k,index2,h,xa,xb,ya,yb,za,zb &
      !                    ,phi,w,x,y,z,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
      !                    ,cocc1,cocc2,cocc3,curv,x1,y1,z1,wnin,kdm,kzone)
      !              wt = (wnout - wnin)*abs(coca3)
      !              a1 = (z1-z(k-1))/hz
      !              b1 = 1.d0 - a1
      !              twz = -((w(i,j,k)-w(i,j,k-1)-tt)/hz+a1*wt)

      !           endif
               else if(phi(i,j,k)*phi(i,j,k+1) <= 0.0) then
                  tuz = -(u(i,j,k)-u(i,j,k-1))/hz
                  tvz = -(v(i,j,k)-v(i,j,k-1))/hz
                  twz = -(w(i,j,k)-w(i,j,k-1))/hz
                  !              elseif(phi(i,j)*phi(i-1,j) .le. 0.0 .or.
                  !    1                phi(i,j)*phi(i+1,j) .le. 0.0) then
               else
                  tuz = -(u(i,j,k+1)-u(i,j,k-1))/hz2
                  tvz = -(v(i,j,k+1)-v(i,j,k-1))/hz2
                  twz = -(w(i,j,k+1)-w(i,j,k-1))/hz2
               end if

               bv(i,j,k) = -( tux*tux + tvy*tvy + twz*twz &
                     + 2.d0*tvx*tuy + 2.d0*twx*tuz + 2.d0*twy*tvz )
               endif !phi(i,j,k) .gt.0  
            end do  !   i= 1,m-1
         end do  !   j= 1,n-1
      end do  !   k= 1,l-1

     if(nse_debug .le.4) then
      do  k= 1,l-1
         do  j= 1,n-1
            do  i= 1,m-1
                
                 
              tux=vey(t,rnu,y(j),x(i),z(k))
              tuy=uey(t,rnu,y(j),x(i),z(k))
              tuz=wey(t,rnu,y(j),x(i),z(k))   
              tvx=uey(t,rnu,x(i),y(j),z(k))
              tvy=vey(t,rnu,x(i),y(j),z(k))
              tvz=wey(t,rnu,x(i),y(j),z(k))   
              twx=uey(t,rnu,x(i),z(k),y(j))
              twy=wey(t,rnu,x(i),z(k),y(j))
              twz=vey(t,rnu,x(i),z(k),y(j))   
               bv(i,j,k) = -( tux*tux + tvy*tvy + twz*twz &
                     + 2.d0*tvx*tuy + 2.d0*twx*tuz + 2.d0*twy*tvz )
             !  write(222222,*) r0,rnu,x(i),y(j),z(k),ui(i,j,k)
            end do  !   i= 1,m-1
         end do  !   j= 1,n-1
      end do  !   k= 1,l-1
      !       et = 0.d0
         cc = 4.0*r0*rnu-sqrt(16.0d0*rnu**2*r0**2-2.0d0*(r0-r0**3)+1.d-16)
         do i0=1,n2
                   ubc2(i0)=-cc*cc/2.0d0/(r0**4)
         enddo
     endif


    !  call iim_poisson_a(m,n,l,n2,xa,xb,ya,yb, &
    !        za,zb,phi,bv,x,y,z,p, elmbda,ubc2,t,rnu)
    !  return
      allocate(fn(n2),fvec(n2),bf(n2),res(n2))
      do i=1,n2
         fn(i) = 0.0d0                   !# start with zero to get -b
         bf(i) = 0.0d0
      end do
 ! cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
 !  do i=1, n2
 !     fn(i)=3*cc*cc/(r0**5)
 !  end do   
         
       ifupdate = .true.
       call matvecp(m,n,l,np,n2,indexx,indexy,indexz &
          ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p,phi,cinfox,cinfoy,cinfoz&
          ,index2,cinfo,fn,fvec,bf,ubc2,bv)
      ifupdate = .false.
      call nsetimer_stop(10)
      call nsetimer_start(8)
      call mpi_barrier(mpi_comm_world,ierr)
      do i=1,n2
         bf(i) = -fvec(i)
      end do
     allocate (desca(9),descb(9))
      ! decompose the matrix to each computer node
    iavg = n2/numtasks
    ires = n2 - iavg*numtasks
    if ( ires == 0 ) then
       ifirst = iavg*mytaskid + 1
       ilast  = iavg*mytaskid + iavg
    else
       if ( (iavg+1)*(mytaskid+1) .le. n2  ) then
          ifirst = (iavg+1)*mytaskid + 1
          ilast  = (iavg+1)*mytaskid + iavg + 1
       else
          ifirst = (iavg+1)*mytaskid+1
          ilast  = n2
       end if
    end if
    nprow=1;npcol=numtasks;
    if(ires ==0) then
    mblock=iavg
    else
    mblock=iavg+1
    endif
    nblock=mblock
     
    call mpi_barrier(mpi_comm_world,ierr)
    call sl_init(ictxtA,nprow,npcol)
    call BLACS_GRIDINFO( ictxtA, nprow, npcol, myrow, mycol )
    nmaxr=numroc(n2,mblock,myrow,0,nprow)
    nmaxc=numroc(n2,nblock,myrow,0,npcol)

    allocate(gmat(nmaxr,nmaxc),ipvt(mblock+nmaxr))
    call descinit( desca,n2,n2, mblock,nblock, 0, 0, ictxtA, max(1,nmaxr),infot )
    call descinit( descb,n2,1, nblock, 1, 0, 0, ictxtA, max(1,nmaxr), infot )
     if(ifirst .le. ilast) then
    do i=ifirst,ilast
     do k=1,n2
        fn(k) = 0.0d0
     end do
     fn(i) = float(10*n2*n2)
   
       call matvecp(m,n,l,np,n2,indexx,indexy,indexz &
          ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p,phi,cinfox,cinfoy,cinfoz&
          ,index2,cinfo,fn,fvec,bf,ubc2,bv)
      do k=1,n2
          gmat(k,i-ifirst+1)=fvec(k)/float(10*n2*n2)
      end do

    end do
    endif
    do i=1,n2
      fn(i) = bf(i)
    end do
      call nsetimer_stop(8)
      call nsetimer_start(9)    
      call mpi_barrier(mpi_comm_world,ierr)
   
      call pdgesv(n2,1,gmat,1,1,desca,ipvt,fn,1,1,descb,infot)
      
      call mpi_barrier(mpi_comm_world,ierr)
    call mpi_bcast(fn,n2,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
   call blacs_gridexit(ictxtA)
   !call blacs_exit(1)
!  cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
!   do i=1, n2
!      fn(i)=3*cc*cc/(r0**5)
!   end do   
!   do i=1,n2
!   write(15,*) dot_product(gmatA(i,1:n2),fn(1:n2)),bf(i)
!   enddo
    deallocate(gmat,ipvt,desca,descb)
       call nsetimer_stop(9)
       call nsetimer_start(10)
       call matvecp(m,n,l,np,n2,indexx,indexy,indexz &
          ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p,phi,cinfox,cinfoy,cinfoz&
          ,index2,cinfo,fn,fvec,bf,ubc2,bv)
      call resid(n2,fvec,bf,res)
      err=dsqrt(prod1(n2,res,res))/dsqrt(prod1(n2,bf,bf))
      if(master) then
      write(6,*) err,dsqrt(prod1(n2,res,res)),dsqrt(prod1(n2,bf,bf))
      endif
      if(err .gt. 1.0d-2) then
      if(master) then
      write(6,*) "There is deviation in pressure direct solver "
      endif
       endif
!       if(master) then
!        do i=1,n2
!        write(12,*) ubc2(i),fn(i),fvec(i),bf(i) 
!        enddo   
!        endif           
    deallocate(fn,fvec,bf,res)
      call nsetimer_stop(10) 
!     do k=1,l-1
!        do j=1,n-1
!           do  i=1,m-1
!              index(i,j,k) = 0     !#  regular grid point
!              !             index2(i,j,k) = 0
!              if(phi(i,j,k) == 0.0) index(i,j,k) = 3
!!              if(phi(i,j,k) > 0.0) then
!                 index(i,j,k)=5
!                 if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 4
!!                 if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 4
!                 if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 4

!              end if

!              if(phi(i,j,k) < 0.0) then
!                 index(i,j,k)=1
!                 if(phi(i,j,k)*phi(i-1,j,k) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i+1,j,k) <= 0.0) index(i,j,k) = 2
!                 if(phi(i,j,k)*phi(i,j-1,k) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i,j+1,k) <= 0.0) index(i,j,k) = 2
!                 if(phi(i,j,k)*phi(i,j,k-1) <= 0.0 .or. &
!                       phi(i,j,k)*phi(i,j,k+1) <= 0.0) index(i,j,k) = 2

!              end if
!           end do
!        end do
!     end do
!     do k=1,l-1
!        do j=1,n-1
!           do i=1,m-1
!              write(123,*) p(i,j,k),pe(rnu,t,x(i),y(j),z(k)),index(i,j,k)
!              write(456,*) bv(i,j,k)
!                        write(123789123,*) index(i,j,k),i,j,k,bv(i,j,k),bv(i,j,k)-(p(i-1,j,k)+p(i+1,j,k)+p(i,j+1,k)&
!                       +p(i,j-1,k)+p(i,j,k+1)+p(i,j,k-1)-6*p(i,j,k))/(h*h),&
!                    1-(p(i-1,j,k)+p(i+1,j,k)+p(i,j+1,k) +p(i,j-1,k)+p(i,j,k+1)+p(i,j,k-1)-6*p(i,j,k))/(h*h*bv(i,j,k))

!           end do
!        end do
!     end do
!     ! explicit A for augmented variables
!     ! compute u,v,w





      return
   end subroutine




   !-----------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine project here]
   subroutine project2(m,n,l,hx,hy,hz,i,j,k,x,y,z,x1,y1,z1,phi)
      implicit none

      !-- subroutine project find the projection of the interface of a given
      !   grid point (x_i,y_j)
      
      !   Input: m,n,h,i,j,x,y,phii
      !           see the description in the subroutine iim_poisson.f
      !   Output: x1,y1 the projection on the interface of (x(i),y(j)).

      _REAL_ phi(0:m,0:n,0:l),x(0:m),y(0:n),z(0:l)
      integer m,n,l,i,j,k,info
      _REAL_ x1,y1,z1,hx,hy,hz,a0,a1,a2,dproj,hmax,ph,phn1,phx,phn
      _REAL_ phxx,phxy,phxz,phy,phyy,phyz,phz,phzz,r1,r2,r3,r4,r5
      



      !   Compute the deepest descent direction (-grad phi(i,j)), then
      !   locates the control points on the level set ----------------------

      hmax = max(hx,hy,hz)
      call phiinf(m,n,l,i,j,k,hx,hy,hz,phi,phx,phy,phz &
            ,phxx,phyy,phzz,phxy,phxz,phyz)

      phn = phx*phx + phy*phy + phz*phz
      phn1 = dsqrt(phn + 1.0d-10)
      a0 = phi(i,j,k)
      a1 = -phn
      a2 =(phx*phx*phxx+phy*phy*phyy+phz*phz*phzz)*0.5d0 &
            +phx*phy*phxy+phx*phz*phxz+phy*phz*phyz
      call rootp22(a0,a1,a2,r1,r2,info)

      !-------------------------------------------------
      if(info > 0) then    !# There is at least one root

         r3 = min(r1,r2)*1.0
         r4 = max(r1,r2)*1.0

         if(phi(i,j,k) >= 0.0) then
            if(r3 >= 0.0) then
               r5 = r3
            else
               r5 = r4
            end if
         else
            if(r4 <= 0.0) then
               r5 = r4
            else
               r5 = r3
            end if
         end if
         x1 = x(i) - r5*phx
         y1 = y(j) - r5*phy
         z1 = z(k) - r5*phz

         dproj = (x1-x(i))*(x1-x(i)) + (y1-y(j))*(y1-y(j)) &
               + (z1-z(k))*(z1-z(k))

      end if

      !------------ if no real root is found, or it is too far away (h),
      !     function may be too bad, choose the intersection between
      !     the interface and the grid lines           -------------------

      if( info == 0 .or. sqrt(dproj) > 2.0*hmax) then
         if(phi(i,j,k)*phi(i+1,j,k) < 0.0) then
            y1 = y(j)
            z1 = z(k)
            x1 = x(i) + hx*phi(i,j,k)/(phi(i,j,k)-phi(i+1,j,k))
            return
         end if
         if(phi(i,j,k)*phi(i-1,j,k) < 0.0) then
            y1 = y(j)
            z1 = z(k)
            x1 = x(i) - hx*phi(i,j,k)/(phi(i,j,k)-phi(i-1,j,k))
            return
         end if

         if(phi(i,j,k)*phi(i,j+1,k) < 0.0) then
            x1 = x(i)
            z1 = z(k)
            y1 = y(j) + hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j+1,k))
            return
         end if
         if(phi(i,j,k)*phi(i,j-1,k) < 0.0) then
            x1 = x(i)
            z1 = z(k)
            y1 = y(j) - hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j-1,k))
            return
         end if

         if(phi(i,j,k)*phi(i,j,k+1) < 0.0) then
            x1 = x(i)
            y1 = y(j)
            z1 = z(k) + hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k+1))
            return
         end if
         if(phi(i,j,k)*phi(i,j,k-1) < 0.0) then
            x1 = x(i)
            y1 = y(j)
            z1 = z(k) - hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k-1))
            return
         end if

      end if  ! ( info == 0 .or. dproj > 2.0*hmax)

      !----------------- end -----------------------------------------

      return
   end subroutine

   ! subroutine to reinitialize the level set function to signed distance function.
   ! nx, ny: number of grid points in x and y direction.
   ! dx, dy: grid size in x and y direction.
   ! l:  the lth level set function for the lth domain.
   ! no: total number of subdomains (level set functions).
   ! u:  the level set functions.
   ! fi: working array for reinitialization.
   ! sign: the flag of the band to be reinitialized.
   ! width: the width of the band to be reinitialized.
   ! num: number of iterations. width=0.1 X num X gridsize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine reinit here]
   subroutine reinit(nx,ny,nz,al,dx,dy,dz,l,no,u,fi,sign,width,num)
      implicit none
       integer nx, ny, nz, l, no, num
      real*8  u(nx+1,ny+1,nz+1,no), fi(2,nx+1,ny+1,nz+1)
      integer sign(nx+1,ny+1,nz+1)
      real*8  dx, dy, dz,width,dlx,drx, d1x, d2x,dly, dry, d1y, d2y
      real*8  al, dis, dlz,drz, d1z, d2z
      integer ip1, im1, jp1, jm1
      integer i,j,k,old,new,p,kk,km1,kp1

      !       al=0.0d0
      old=1
      new=2

      ! initialize the data{
      do k=1,nz+1
         if(k == 1) then
            km1=k
         else
            km1=k-1
         end if
         if(k == nz+1) then
            kp1=k
         else
            kp1=k+1
         end if
         do j=1,ny+1
            if(j == 1) then
               jm1=j
            else
               jm1=j-1
            end if
            if(j == ny+1) then
               jp1=j
            else
               jp1=j+1
            end if
            do i=1,nx+1
               if(i == 1) then
                  im1=i
               else
                  im1=i-1
               end if
               if(i == nx+1) then
                  ip1=i
               else
                  ip1=i+1
               end if
               dis=min(dabs(u(im1,j,k,l)), dabs(u(ip1,j,k,l)), &
                     dabs(u(i,jm1,k,l)), dabs(u(i,jp1,k,l)), &
                     dabs(u(i,j,km1,l)), dabs(u(i,j,kp1,l)), &
                     dabs(u(i,j,k,l)))
               if(dis < width) then
                  sign(i,j,k)=1
               else
                  sign(i,j,k)=0
               end if
               if(dabs(u(i,j,k,l)) > width) then
                  if(u(i,j,k,l) > 0.0d0) then
                     fi(old,i,j,k)=width
                  else
                     fi(old,i,j,k)=-width
                  end if
               else
                  fi(old,i,j,k)=u(i,j,k,l)
               end if
            end do
         end do  !  j=1,ny+1
      end do  !  k=1,nz+1
      ! end if initialization}
      
      ! do iterations to reinitialize{
      do kk=1,num
         do k=2,nz
            do j=2,ny
               do i=2,nx
                  if(sign(i,j,k) == 1) then
                     ! compute the derivative in x using 2nd ENO{
                     if (i == 2) then
                        call eno2(0.0d0,fi(old,i-1,j,k),fi(old,i,j,k), &
                              fi(old,i+1,j,k),fi(old,i+2,j,k), dx, i,nx, &
                              dlx,drx, d1x, d2x)
                     end if
                     if (i > 2 .and. i < nx) then
                        call eno2(fi(old,i-2,j,k),fi(old,i-1,j,k),fi(old,i,j,k), &
                              fi(old,i+1,j,k),fi(old,i+2,j,k), dx, i,nx, &
                              dlx,drx, d1x, d2x)
                     end if
                     if (i == nx) then
                        call eno2(fi(old,i-2,j,k),fi(old,i-1,j,k),fi(old,i,j,k), &
                              fi(old,i+1,j,k), 0.0d0,     dx, i,nx, &
                              dlx,drx, d1x, d2x)
                     end if
                     ! end of computing the derivative in x using 2nd ENO}

                     ! compute the derivative in y using 2nd ENO{
                     if (j == 2) then
                        call eno2( 0.0d0, fi(old,i,j-1,k), &
                              fi(old,i,j,k), fi(old,i,j+1,k), fi(old,i,j+2,k), &
                              dy, j, ny, dly, dry, d1y, d2y)
                     end if
                     if (j > 2 .and. j < ny) then
                        call eno2( fi(old,i,j-2,k), fi(old,i,j-1,k), &
                              fi(old,i,j,k), fi(old,i,j+1,k), fi(old,i,j+2,k), &
                              dy, j, ny, dly, dry, d1y, d2y)
                     end if
                     if (j == ny) then
                        call eno2( fi(old,i,j-2,k), fi(old,i,j-1,k), &
                              fi(old,i,j,k), fi(old,i,j+1,k), 0.0d0, &
                              dy, j, ny, dly, dry, d1y, d2y)
                     end if
                     ! end of computing the derivative in y using 2nd ENO}

                     ! compute the derivative in z using 2nd ENO{
                     if (k == 2) then
                        call eno2( 0.0d0, fi(old,i,j,k-1), &
                              fi(old,i,j,k), fi(old,i,j,k+1), fi(old,i,j,k+2), &
                              dz, k, nz, dlz, drz, d1z, d2z)
                     end if
                     if (k > 2 .and. k < nz) then
                        call eno2( fi(old,i,j,k-2), fi(old,i,j,k-1), &
                              fi(old,i,j,k), fi(old,i,j,k+1), fi(old,i,j,k+2), &
                              dz, k, nz, dlz, drz, d1z, d2z)
                     end if
                     if (k == nz) then
                        call eno2( fi(old,i,j,k-2), fi(old,i,j,k-1), &
                              fi(old,i,j,k), fi(old,i,j,k+1), 0.0d0, &
                              dz, k, nz, dlz, drz, d1z, d2z)
                     end if
                     ! end of computing the derivative in y using 2nd ENO}

                     ! update the distance function{
                     if(u(i,j,k,l) < -al) then
                        fi(new,i,j,k)=fi(old,i,j,k)-0.1d0*dx* &
                              (u(i,j,k,l)/dsqrt(u(i,j,k,l)**2+dx))* &
                              (dsqrt((min(dlx,0.0d0))**2+(max(drx,0.0d0))**2+ &
                              (min(dly,0.0d0))**2+(max(dry,0.0d0))**2+ &
                              (min(dlz,0.0d0))**2+(max(drz,0.0d0))**2) &
                              -1.0d0)
                     else if(u(i,j,k,l) > al) then
                        fi(new,i,j,k)=fi(old,i,j,k)-0.1d0*dx* &
                              (u(i,j,k,l)/dsqrt(u(i,j,k,l)**2+dx))* &
                              (dsqrt((max(dlx,0.0d0))**2+(min(drx,0.0d0))**2+ &
                              (max(dly,0.0d0))**2+(min(dry,0.0d0))**2+ &
                              (max(dlz,0.0d0))**2+(min(drz,0.0d0))**2) &
                              -1.0d0)
                     else
                        fi(new,i,j,k)=fi(old,i,j,k)
                     end if

                  else
                     fi(new,i,j,k)=fi(old,i,j,k)
                  end if  ! (sign(i,j,k) == 1)
                  ! end of update the distance function}
               end do  !  i=2,nx
            end do  !  j=2,ny
         end do  !  k=2,nz

         ! boundary condition{
         ! Set extropolation boundary condition
         call ex_bo1(fi, 1, nx, ny, nz, old, new)
         ! end of boundary condition}
         p=new
         new=old
         old=p
      end do  !  kk=1,num
      ! end of reinitialization}

      do k=1,nz+1
         do j=1,ny+1
            do i=1,nx+1
               if(sign(i,j,k) == 1) then
                  u(i,j,k,l)=fi(old,i,j,k)
               end if
            end do
         end do
      end do

      return
   end subroutine

   !----------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine eno2 here]
   subroutine eno2(um2,um1,ui,up1,up2, h, i,n, &
         dum,dup, duc, du2)
      integer i, n
      real*8 um2,um1,ui,up1,up2,    h
      real*8 dum, dup, duc, du2
      !       real*8 absmin

      real*8 dl1, dl2, dr1, dr2

      if(i > 2) then
         dl1=(ui  - 2.0d0*um1 + um2)/h**2
         dl2=(up1 - 2.0d0*ui  + um1)/h**2
      else
         dl1=0.0d0
         dl2=0.0d0
      end if
      dum=(ui - um1)/h+0.5d0*h*absmin(dl1,dl2)

      if(i < n) then
         dr1=(up2 - 2.0d0*up1 +ui)/h**2
         dr2=(up1 - 2.0d0*ui + um1)/h**2
      else
         dr1=0.0d0
         dr2=0.0d0
      end if
      dup=(up1 - ui)/h-0.5d0*h*absmin(dr1, dr2)

      duc=0.5d0*(up1 - um1)/h
      du2=(up1 - 2.0d0*ui + um1)/h**2

      return
   end subroutine
   !--------------------------------------

   !   Set boundary by extropolation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ex_bo1 here]
   subroutine ex_bo1(u, no, nx, ny, nz, old, new)
      !       -----------------------------------------
      integer no, old, new, nx, ny, nz
      real*8 u(no, 2, nx+1, ny+1, nz+1)

      integer l, i, j, k
      ! boundary
      do l=1,no
         do j=2,ny
            do i=2,nx
               u(l,new,i,j,1)=2.0d0*u(l,new,i,j,2)-u(l,new,i,j,3)
               u(l,new,i,j,nz+1)=2.0d0*u(l,new,i,j,nz)-u(l,new,i,j,nz-1)
            end do
         end do
         do k=2,nz
            do i=2,nx
               u(l,new,i,1,k)=2.0d0*u(l,new,i,2,k)-u(l,new,i,3,k)
               u(l,new,i,ny+1,k)=2.0d0*u(l,new,i,ny,k)-u(l,new,i,ny-1,k)
            end do
         end do
         do k=2,nz
            do j=2,ny
               u(l,new,1,j,k)=2.0d0*u(l,new,2,j,k)-u(l,new,3,j,k)
               u(l,new,nx+1,j,k)=2.0d0*u(l,new,nx,j,k)-u(l,new,nx-1,j,k)
            end do
         end do

         do i=2,nx
            u(l,new,i,1,1)=2.0d0*u(l,new,i,2,2)-u(l,new,i,3,3)
            u(l,new,i,ny+1,1)=2.0d0*u(l,new,i,ny,2)-u(l,new,i,ny-1,3)
            u(l,new,i,1,nz+1)=2.0d0*u(l,new,i,2,nz)-u(l,new,i,3,nz-1)
            u(l,new,i,ny+1,nz+1)=2.0d0*u(l,new,i,ny,nz) &
                  -u(l,new,i,ny-1,nz-1)
         end do

         do j=2,ny
            u(l,new,1,j,1)=2.0d0*u(l,new,2,j,2)-u(l,new,3,j,3)
            u(l,new,nx+1,j,1)=2.0d0*u(l,new,nx,j,2)-u(l,new,nx-1,j,3)
            u(l,new,1,j,nz+1)=2.0d0*u(l,new,2,j,nz)-u(l,new,3,j,nz-1)
            u(l,new,nx+1,j,nz+1)=2.0d0*u(l,new,nx,j,nz) &
                  -u(l,new,nx-1,j,nz-1)
         end do

         do k=2,nz
            u(l,new,1,1,k)=2.0d0*u(l,new,2,2,k)-u(l,new,3,3,k)
            u(l,new,nx+1,1,k)=2.0d0*u(l,new,nx,2,k)-u(l,new,nx-1,3,k)
            u(l,new,1,ny+1,k)=2.0d0*u(l,new,2,ny,k)-u(l,new,3,ny-1,k)
            u(l,new,nx+1,ny+1,k)=2.0d0*u(l,new,ny,nz,k) &
                  -u(l,new,ny-1,nz-1,k)
         end do

         u(l,new,1,1,1)=2.0d0*u(l,new,2,2,2)-u(l,new,3,3,3)
         u(l,new,nx+1,1,1)=2.0d0*u(l,new,nx,2,2)-u(l,new,nx-1,3,3)
         u(l,new,1,ny+1,1)=2.0d0*u(l,new,2,ny,2)-u(l,new,3,ny-1,3)
         u(l,new,nx+1,ny+1,1)=2.0d0*u(l,new,nx,ny,2) &
               -u(l,new,nx-1,ny-1,3)
         u(l,new,1,1,nz+1)=2.0d0*u(l,new,2,2,nz)-u(l,new,3,3,nz-1)
         u(l,new,nx+1,1,nz+1)=2.0d0*u(l,new,nx,2,nz) &
               -u(l,new,nx-1,3,nz-1)
         u(l,new,1,ny+1,nz+1)=2.0d0*u(l,new,2,ny,nz) &
               -u(l,new,3,ny-1,nz-1)
         u(l,new,nx+1,ny+1,nz+1)=2.0d0*u(l,new,nx,ny,nz) &
               -u(l,new,nx-1,ny-1,nz-1)
      end do  !  l=1,no

      return
   end subroutine

   !---------------------------------------------------------------

   real*8 function sgn(x)
   real*8 x

   if (x >= 0.0) then
      sgn=1.0d0
   else
      sgn=-1.0d0
   end if

   return
end function
!FIXME:                nse.f, line 8163: bad indentation level at end of subroutine ex_bo1
   !---------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine weno3 here]
   subroutine weno3(um2, um1, ui, up1, up2, h, i, n, &
         dum, dup, duc, du2)

      integer i, n
      real*8  um2,um1,ui,up1,up2, h
      real*8 dum, dup, duc, du2
      real*8 wm, wp
      real*8 d1, d2, d3, d4, dd1, dd2, dd

      real*8 ep, r

      ep=1.0d-6
      if(i == 2) then
         um2=um1+sgn(um1)*dabs(um1-ui)
      else if(i == n) then
         up2=up1+sgn(up1)*dabs(up1-ui)
      end if

      d1=(um1-um2)/h
      d2=(ui-um1)/h
      d3=(up1-ui)/h
      d4=(up2-up1)/h
      dd1=um2-2.0d0*um1+ui
      dd=um1-2.0d0*ui+up1
      dd2=ui-2.0d0*up1+up2

      r=(ep+dd1*dd1)/(ep+dd*dd)
      wm=1.0d0/(1.0d0+2.0d0*r*r)
      r=(ep+dd2*dd2)/(ep+dd*dd)
      wp=1.0d0/(1.0d0+2.0d0*r*r)

      duc=0.5d0*(d2+d3)
      dum=duc-0.5d0*wm*(d1-2.0d0*d2+d3)
      dup=duc-0.5d0*wp*(d2-2.0d0*d3+d4)
      du2=(d3-d2)/h

      return
   end subroutine

   !---------------------------------------------------------------



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rootp3 here]
   subroutine rootp3(a0,a1,a2,a3,r1,r2,r3,tol,info)

      !*******************************************************************************
      
      !  rootp3 finds the real root(s) for cubic polynomial
      !      a3 x^3 + a2 x^2 + a1 x + a0
      !  in double precision. This subroutine uses Newton's method to find
      !  the maximum root, then use deflation technique to get the quadratic
      !  polynomial. At last this routine uses rootp2 to find other roots.
      
      !  a0,a1,a2,a3: On entry, the coefficients of the cubic.
      !  info: Output which indicates how many real roots are found. info takes
      !        value 0,1,2,3, or 100. info = 100 means a3=a2=a1=0, so the input
      !        coefficients are wrong. info = 0 means a3 = 0, so the cubic is
      !        actually a quadratic and no real root is found. info = 2 means a3 = 0,
      !        but there are two real roots.
      !  r1,r2,r3, Output: The real roots found. some of r3,r2,r1 can be dumb if
      !        info < 3.
      
      !*******************************************************************************

      implicit none
      _REAL_ a0,a1,a2,a3,r1,r2,r3,r4,rm,tol,a,b
      integer info,k

      if( abs(a3) <= 1.0d-10) then
         call rootp22(a0,a1,a2,r1,r2,info)
         return
      else
         r4 = max( abs(a0/a3), abs(a1/a3)+1.0d0, abs(a2/a3)+1.0d0)
         a = r4
         b = -r4

         !------ Bisection to get good initial guess ----------------------------------

         do k=1,10
            rm = 0.5d0*(a+b)
            if( f(a0,a1,a2,a3,rm)*f(a0,a1,a2,a3,a) <= 0.0d0) then
               b = rm
            else
               a = rm
            end if
         end do

         r4 = 0.5d0*(a+b)

         10 r1 = r4 - f(a0,a1,a2,a3,r4)/d(a0,a1,a2,a3,r4)
         if( abs(r4-r1) > tol .and. &
               abs(f(a0,a1,a2,a3,r4)) > tol) then
            r4 = r1
            goto 10
         else
            a2 = a3*r1 + a2
            a1 = a2*r1 + a1
            call rootp22(a1,a2,a3,r2,r3,info)
            info = info + 1
            return
         end if
      end if

   end subroutine

   !----------------------------------------------------------------------------

   double precision function f(a0,a1,a2,a3,x)
   implicit none
   _REAL_ f,x,a0,a1,a2,a3
   f = a0 + x*(a1+x*(a2+a3*x))
   return
end function
!FIXME:                nse.f, line 8276: bad indentation level at end of subroutine rootp3

   double precision function d(a0,a1,a2,a3,x)
   implicit none
   _REAL_ a0,a1,a2,a3,x  
   d = a1 + x*(2.0d0*a2 +x*3.0d0*a3)
   return
end function
!FIXME:                nse.f, line 8283: bad indentation level at end of subroutine rootp3

   !------------------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rootp2 here]
   subroutine rootp22(a0,a1,a2,r1,r2,info)
      implicit none
      _REAL_ a0,a1,a2,r1,r2,t
      integer info
       

      if( abs(a2) <= 1.0d-10) then
         if( abs(a1) <= 1.0d-10) then
            if( abs(a0) <= 1.0d-10) then
               info = 10
            else
               info = 100
            end if
            return
         else
            info = 2
            r1 = -a0/a1
            r2 = r1
            return
         end if
      else
         t = a1*a1 - 4.0d0*a0*a2
         if ( t >= 0.0) then
            info = 2
            t = dsqrt(t)
            if( a1 >= 0.0d0) then
               r1 = ( -a1 - t)/(2.0d0*a2)
            else
               r1 = ( -a1 + t)/(2.0d0*a2)
            end if

            if( t == 0.0) then
               r2 = r1
            else
               r2 = a0/(a2*r1)
            end if
            return
         else
            info = 0
            return
         end if
      end if  ! ( abs(a2) <= 1.0d-10)

   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function ue here]
   function ue(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,ue,cc,r,rs
!      ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)

      if(rs <= r0) then
         ue = cc*x/(r0**3)
         !            ue=0
      else
         ue = cc*x/(rs*r)
      end if

      ue = ue*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function ve here]
   function ve(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,ve,cc,r,rs
  !    ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(rs <= r0) then
         ve = cc*y/(r0**3)
         !           ve=0
      else
         ve = cc*y/(rs*r)
      end if

      ve = ve*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function we here]
   function we(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,we,cc,r,rs
   !   ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(rs <= r0) then
         we = cc*z/(r0**3)
         !           we=0
      else
         we = cc*z/(rs*r)
      end if

      we = we*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function uey here]
   function uey(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,uey,cc,r,rs
   !   ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)

      if(rs <= r0) then
         uey = 0.0d0
      else
         uey = -3.0d0*cc*x*y/(r*r*rs)
      end if

      uey = uey*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function vey here]
   function vey(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,vey,cc,r,rs
   !   ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(rs <= r0) then
         vey = cc/(r0**3)
      else
         vey = cc*(x*x+z*z-2.0*y*y)/(r*r*rs)
      end if

      vey = vey*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function wey here]
   function wey(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,wey,cc,r,rs
  !    ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(rs <= r0) then
         wey = 0.d0
      else
         wey = -3.d0*cc*z*y/(r*r*rs)
      end if

      wey = wey*fti(t)

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function fe here]
   function fe(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,fe,r,rs,pi

   !   ! common /radius/r0

      pi = 4.*datan(1.0d0)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(r < r0*r0) then
         !          uet = y*(r - 0.25)
         !          uconv = 2.0*y*y*(r - 0.25)*x -
         !    1              x*(r - 0.25)*(x*x+3.0*y*y-0.25)
         !          ulap = 8.0*rnu*y                  !# Laplacian term
         !          fe = uet*dfti(t) + fti(t)*fti(t)*uconv - fti(t)*ulap
         fe = 0.0
      else
         !          fet = (-1.0+2.0*rs)
         !          uconv = -x*fet*fet/r
         !          ulap = - rnu*y/(rs*r)
         !          uet = y/rs - 2.0*y
         !          fe =uet*dfti(t)+fti(t)*fti(t)*uconv-fti(t)*ulap
         !          fe = fe -pi*dsin(pi*x)*dcos(pi*y)*0.0
         fe = 0.0
      end if

      return
   end function

   !# *******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function ge here]
   function ge(t,rnu,x,y,z)
      implicit none
      _REAL_ rnu,t,x,y,z,ge,r,rs,pi
   !   ! common /radius/r0

      pi = 4.*datan(1.0d0)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(r < r0*r0 ) then
         !          vet = -x*(r-0.25)
         !          vconv = y*(r-0.25)*(-3.0*x*x-y*y+0.25)+
         !    1         2.0*x*x*(r-0.25)*y
         !          vlap = -8.0*rnu*x
         !          ge = vet*dfti(t)+fti(t)*fti(t)*vconv-fti(t)*vlap
         ge = 0.0
      else
         !          get = (-1.0+2.0*rs)
         !          vconv = -y*get*get/r
         !          vlap = rnu*x/(rs*r)
         !          vet = -x/rs+ 2.0*x
         !          ge = vet*dfti(t)+fti(t)*fti(t)*vconv-fti(t)*vlap
         !          ge = ge -pi*dsin(pi*y)*dcos(pi*x)*0.0
         ge = 0.0

      end if

      return
   end function

   !# *********************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function he here]
   function he(t,rnu,x,y,z)
      implicit none
      
      _REAL_ rnu,t,x,y,z,he,pi,r,rs

   !   ! common /radius/r0

      pi = 4.*datan(1.0d0)
      r = x*x + y*y + z*z
      rs = dsqrt(r)
      if(r < r0*r0) then
         he = 0.0
      else
         he = 0.0
      end if

      return
   end function

   !# ************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function fti here]
   function fti(t)
      implicit none
      _REAL_ fti,t

      !       fti = dsin(t)
      fti = 1.0
      !       fti = dexp(-t)
      !       fti = 1.0 + dexp(-t)

      return
   end function

   !# ******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function dfti here]
   function dfti(t)
      implicit none
      _REAL_ dfti,t

      !       dfti = dcos(t)
      dfti = 0.0
      !       dfti = -dexp(-t)

      return
   end function

   !# *******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function pe here]
   function pe(rnu,t,x,y,z)
      implicit none
      _REAL_ pe,rnu,t,x,y,z,cc,pi,r,r2,rs
   !   ! common /radius/r0

      cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      pi = 4.*datan(1.0d0)

      r = x*x + y*y + z*z
      rs = sqrt(r)
      r2 = r0*r0

      if(rs <= r0) then
         pe = -0.5*cc*cc*r/(r2**3)
         !          pe=0
      else
         pe = -0.5*cc*cc/(r*r)
      end if


      return
   end function

   !# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiset here]
   subroutine phiset(m,n,l,r0,x0,y0,z0,x,y,z,phi)
      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      _REAL_ r0,a0,b0,c0,xt,yt,zt,x0,y0,z0
      integer m,n,l,i,j,k
 !      ! common /radius/r0

      a0 = r0
      b0 = r0
      c0 = r0
      do k=0,l
         do j=0,n
            do i=0,m

               xt = x(i) + 0.0d0-x0
               yt = y(j) + 0.0d0-y0
               zt = z(k) + 0.0d0-z0
               phi(i,j,k) = dsqrt(xt*xt+yt*yt &
                     +zt*zt)-r0

            end do
         end do
      end do

      return
   end subroutine

   !# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiset2 here]
   subroutine phiset3(m,n,l,x,y,z,phi)
      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      _REAL_ r0,a0,b0,c0,xt,yt,zt,xt1,xt2,xt3,yt1,yt2,yt3,zt1,zt2,zt3,phi1,phi2,phi3
      integer m,n,l,i,j,k
!      ! common /radius/r0

      a0 = r0
      b0 = r0
      c0 = r0

      !       print *,'r0=',r0
      zonexyz(1,1)=0.732d0  
      zonexyz(2,1)=-0.719d0
      zonexyz(3,1)=-0.717d0
      zonexyz(1,2)=0.743d0
      zonexyz(2,2)=0.722d0
      zonexyz(3,2)=-0.734d0
      zonexyz(1,3)=0.051d0
      zonexyz(2,3)=-0.012d0
      zonexyz(3,3)=0.023d0

      do k=0,l
         do j=0,n
            do i=0,m

               !           phi1 = dsqrt( (x(i)+1.2)**2 + y(j)**2 )-r0
               xt = x(i)
               yt = y(j)
               zt = z(k)
               xt1 = x(i)-zonexyz(1,1)
               xt2 = x(i)-zonexyz(2,1)
               xt3= x(i)-zonexyz(3,1) 
               yt1 = y(j)-zonexyz(1,2)
               yt2 = y(j)-zonexyz(2,2)
               yt3= y(j)-zonexyz(3,2)
               zt1= z(k)-zonexyz(1,3)
               zt2 =z(k)-zonexyz(2,3)
               zt3 =z(k)-zonexyz(3,3)

               !           call rose(7,xt,yt,zt,r0,re,phi3)
               call rose(6,xt1,yt1,zt1,0.0d0,phi1)
               call rose(5,xt2,yt2,zt2,0.0d0,phi2)
               call rose(5,xt3,yt3,zt3,0.0d0,phi3)
               !           call rose(6,xt,yt1,r0,0.0d0,phi1)
               !           call rose(5,xt,yt2,re,0.0d0,phi2)
               !           call rose(6,xt,yt1,r0,0.09d0,phi4)
               !            call rose(5,xt,yt2,r0,0.04d0,phi5)
               !           call rose(5,xt,yt2,r0,0.08d0,phi5)

               !        write(*,*)dsqrt( (x(i))**2 + y(j)**2 )-r0,phi3
               !        phi3 = dsqrt( (x(i))**2 + y(j)**2 )-r0

               !        phi2 = dsqrt( (x(i)-1.2)**2 + y(j)**2 )-r0
               !        phi3 = dsqrt( (x(i))**2 + y(j)**2 )-r0
               !        phi4 = dsqrt( (x(i))**2 + (y(j)-1.2)**2 )-r0
               !        phi5 = dsqrt( (x(i))**2 + (y(j)+1.2)**2 )-r0

               !       phi(i,j) = min(phi1,phi2,phi3,phi4,phi5)
               !       phi(i,j) = phi3
               phi(i,j,k) = min(phi1,phi2,phi3)

            end do
         end do
      end do

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phiset2 here]
   subroutine phiset2(m,n,l,x,y,z,phi)
      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      _REAL_ r0,a0,b0,c0,xt,yt,zt,xt1,xt2,yt1,yt2,zt1,zt2,phi1,phi2
      integer m,n,l,i,j,k
!      ! common /radius/r0

      a0 = r0
      b0 = r0
      c0 = r0

      !       print *,'r0=',r0
      

      do k=0,l
         do j=0,n
            do i=0,m

               !           phi1 = dsqrt( (x(i)+1.2)**2 + y(j)**2 )-r0
               xt = x(i)
               yt = y(j)
               zt = z(k)
               xt1 = x(i)+0.71d0
               xt2 = x(i)-0.72d0
               yt1 = y(j)+0.07d0
               yt2 = y(j)-0.10d0
               zt1= z(k)-0.02d0
               zt2 = z(k)+0.05d0 

               !           call rose(7,xt,yt,zt,r0,re,phi3)
               call rose(6,xt1,yt1,zt1,0.0d0,phi1)
               call rose(5,xt2,yt2,zt2,0.0d0,phi2)
               !           call rose(6,xt,yt1,r0,0.0d0,phi1)
               !           call rose(5,xt,yt2,re,0.0d0,phi2)
               !           call rose(6,xt,yt1,r0,0.09d0,phi4)
               !            call rose(5,xt,yt2,r0,0.04d0,phi5)
               !           call rose(5,xt,yt2,r0,0.08d0,phi5)

               !        write(*,*)dsqrt( (x(i))**2 + y(j)**2 )-r0,phi3
               !        phi3 = dsqrt( (x(i))**2 + y(j)**2 )-r0

               !        phi2 = dsqrt( (x(i)-1.2)**2 + y(j)**2 )-r0
               !        phi3 = dsqrt( (x(i))**2 + y(j)**2 )-r0
               !        phi4 = dsqrt( (x(i))**2 + (y(j)-1.2)**2 )-r0
               !        phi5 = dsqrt( (x(i))**2 + (y(j)+1.2)**2 )-r0

               !       phi(i,j) = min(phi1,phi2,phi3,phi4,phi5)
               !       phi(i,j) = phi3
               phi(i,j,k) = min(phi1,phi2)

            end do
         end do
      end do

      return
   end subroutine

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fuj here]
   subroutine fuj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3 &
         ,curv,fujv)
      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      _REAL_ t,rnu,x1,y1,z1,coc1,coc2,coc3,curv,fujv,pi,r,rs
      integer m,n,l         



      pi = 4.*datan(1.0d0)
      r = x1*x1 + y1*y1 + z1*z1
      rs = dsqrt(r)

      !          fet = (-1.0+2.0*rs)
      !          uconv = -x1*fet*fet/r
      !          ulap = - rnu*y1/(rs*r)
      !          uet = y1/rs - 2.0*y1
      !          fujv =uet*dfti(t)+fti(t)*fti(t)*uconv-fti(t)*ulap
      !          fujv =  fujv - pi*dsin(pi*x1)*dcos(pi*y1)*0.0

      fujv = 0.d0

      return
   end subroutine

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fvj here]
   subroutine fvj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3 &
         ,curv,fvjv)
      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      integer m,n,l
      _REAL_ x1,y1,z1,coc1,coc2,coc3,curv,fvjv,t,rnu,pi,r,rs     

      pi = 4.*datan(1.0d0)
      r = x1*x1 + y1*y1 + z1*z1
      rs = dsqrt(r)

      !          get = (-1.0+2.0*rs)
      !          vconv = -y1*get*get/r
      !          vlap = rnu*x1/(rs*r)
      !          vet = -x1/rs+ 2.0*x1
      !          fvjv = vet*dfti(t)+fti(t)*fti(t)*vconv-fti(t)*vlap
      !          fvjv = fvjv - pi*dsin(pi*y1)*dcos(pi*x1)*0.0

      fvjv = 0.d0

      return
   end subroutine

   !#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fwj here]
   subroutine fwj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3 &
         ,curv,fwjv)
      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l)
      integer m,n,l
      _REAL_ t,rnu,coc1,coc2,coc3,curv,fwjv,x1,y1,z1,pi,r,rs     


      pi = 4.*datan(1.0d0)
      r = x1*x1 + y1*y1 + z1*z1
      rs = dsqrt(r)

      !          get = (-1.0+2.0*rs)
      !          vconv = -y1*get*get/r
      !          vlap = rnu*x1/(rs*r)
      !          vet = -x1/rs+ 2.0*x1
      !          fvjv = vet*dfti(t)+fti(t)*fti(t)*vconv-fti(t)*vlap
      !          fvjv = fvjv - pi*dsin(pi*y1)*dcos(pi*x1)*0.0

      fwjv = 0.d0

      return
   end subroutine

   !#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rose here]
   subroutine rose(k,x,y,z,re,phi)
      implicit none
       _REAL_ x,y,z,re,phi,pi
       integer k 

      pi = 3.14159265358979d0

      !           if(abs(x) .lt. 1.0e-13) then
      !              if(y .ge. 0.0d0) then
      !                 theta = pi/2.0
      !              else
      !                 theta = -pi/2.0
      !              endif
      !           else
      !              theta = datan(y/x)
      !              if(x .le. 0.0d0) theta = pi + theta
      !           endif

      !           phi =dsqrt(x*x+y*y) - (r0 + re*dsin(k*theta) )
      phi =dsqrt(x*x+y*y+z*z) - r0

      return
   end subroutine

   ! this function measures the volume of the level set function
   ! u(i,j) is the level set function
   ! volume=\int H(u(i,j))dxdy, H is the numerical Heaviside function.
   ! al is the numerical width of the Heaviside function

   real*8  function vol(u,kzone,idx,al,dx,dy,dz,nx,ny,nz)
   _REAL_  u(nx+1,ny+1,nz+1), dx, dy, dz, al
   _REAL_  da, pi
   integer idx,jdx,kk,km1,kp1
   integer kzone(nx+1,ny+1,nz+1), nzone
   integer nx, ny, nz
   integer i,j,k

   pi=2.0d0*dasin(1.0d0)
   vol=0.0d0
   da=dx*dy*dz
   jdx = idx*100000

   do k=1, nz+1
      do j=1, ny+1
         do i=1, nx+1
            kk = kzone(i,j,k)
            if ( kk == jdx .or. match(kk,idx) == 1 ) &
                  vol=vol+heaviside(-u(i,j,k), al, pi, 1)*da
         end do
      end do
   end do

   return
end function
!FIXME:                nse.f, line 8779: bad indentation level at end of subroutine rose
   !************************************************
   real*8 function heaviside(x,al,pi,order)
   !       ----------------------------------------
   !       this function evaluate numerically the Heaviside function
   !       (integration) of the corresponding delta function
   !       x      is the independent variable
   !       al     the width of the delta function
   !       order  is the corresponding order of the delta function
   
   real*8 x, al, pi
   integer order

   heaviside = 0.0d0

   ! 1st order Heaviside function
   if(order == 1) then
      if(x < -al) then
         heaviside=0.0d0
      else if(dabs(x) <= al) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/2.0d0
      else
         heaviside=1.0d0
      end if
   end if

   ! 2nd order Heaviside function
   if (order == 2) then
      if(x < -al) then
         heaviside=0.0d0
      else if(x >= -al.and.x < -0.5d0*al) then
         heaviside=-(1.0d0+x/al+dsin(pi*x/al)/pi)/6.0d0
      else if(dabs(x) <= 0.5d0*al) then
         heaviside=-(1.0d0+x/al+dsin(pi*x/al)/pi)/6.0d0+ &
               (2.0d0+4.0d0*x/al+2.0d0*dsin(2.0d0*pi*x)/pi)/3.0d0
      else if(x > 0.5d0*al.and.x <= al) then
         heaviside=-(1.0d0+x/al+dsin(pi*x/al)/pi)/6.0d0+4.0d0/3.0d0
      else
         heaviside=1.0d0
      end if
   end if

   ! 3rd order Heaviside function
   if(order == 3) then
      if(x < -al) then
         heaviside=0.0d0
      else if(x >= -al.and.x < -0.5d0*al) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/48.0d0
      else if(x >= -0.5d0*al.and.x < -al/3.0d0) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/48.0d0- &
               4.0d0*(2.0d0+4.0d0*x/al+2.0d0*dsin(2.0d0*pi*x)/pi)/15.0d0
      else if(dabs(x) <= al/3.0d0) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/48.0d0- &
               4.0d0*(2.0d0+4.0d0*x/al+2.0d0*dsin(2.0d0*pi*x)/pi)/15.0d0+ &
               81.0d0*(1.0d0+3.0d0*x/al+dsin(3.0d0*pi*x/al)/pi)/80.0d0
      else if(x > al/3.0d0.and.x <= 0.5d0*al) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/48.0d0- &
               4.0d0*(2.0d0+4.0d0*x/al+2.0d0*dsin(2.0d0*pi*x)/pi)/15.0d0+ &
               81.0d0/40.0d0
      else if(x > 0.5*al.and.x <= al) then
         heaviside=(1.0d0+x/al+dsin(pi*x/al)/pi)/48.0d0- &
               16.0d0/15.0d0+81.0d0/40.0d0
      else
         heaviside=1.0d0
      end if
   end if

   if ( (order > 3) .or. (order < 1)) then
      stop 'wrong order'
   end if

   return
end function
!FIXME:                nse.f, line 8850: bad indentation level at end of subroutine rose



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine xcor here]
   subroutine xcor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z &
         ,xg,x1,y1,z1,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor &
         ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,fuju,fvjv,fwjw)


      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l)
      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,x1,y1,z1,rnu
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,unbj,vnbj,wnbj,uncj,vncj,wncj
      _REAL_ ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ucor,vcor,wcor,unji
      _REAL_ uxj,uxxj,vnji,vxj,vxxj,wnji,wxj,wxxj,fuju,fvjv,fwjw
 
      !       call fuj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fuju)
      !       call fvj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fvjv)
      !       call fwj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fwjw)
      ! above outputs fuju fvjv fwjw all zero --CQ

      unji = ftv
      vnji = fnv
      wnji = fmv

      uxj = unji*coca1
      vxj = vnji*coca1
      wxj = wnji*coca1
      !       uyj = unji*coca2
      !       vyj = vnji*coca2
      !       wyj = wnji*coca2
      !       uzj = unji*coca3
      !       vzj = vnji*coca3
      !       wzj = wnji*coca3

      !       uyyjt = -unji*curv    !# = (ft * tau *curv )/rnu
      !       vyyjt = -vnji*curv
      !       wyyjt = -wnji*curv

      !       uxyjt = dft
      !       vxyjt = dfn
      !       wxyjt = dfm

      !       uxzjt = dft
      !       vxzjt = dfn
      !       wxzjt = dfm

      !       uyzjt = dft
      !       vyzjt = dfn
      !       wyzjt = dfm

      !       call bilinear_m(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1
      !    1    ,x,y,z,phi,u,v,w,uxj,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj
      !    2    ,u1p,v1p,w1p)

      !       uxxjt = -uyyjt + ( - fuju)/rnu
      !       vxxjt = -vyyjt + ( - fvjv)/rnu
      !       wxxjt = -wyyjt + ( - fwjw)/rnu

      !       uzzjt = -uyyjt + ( - fuju)/rnu
      !       vzzjt = -vyyjt + ( - fvjv)/rnu
      !       wzzjt = -wyyjt + ( - fwjw)/rnu

      !       uxxj = uxxjt*coca1*coca1+uyyjt*coca2*coca2+uzzjt*coca3*coca3
      !    1  - 2.0*(uxyjt*coca1*coca2+uxzjt*coca1*coca3+uyzjt*coca2*coca3)
      !       vxxj = vxxjt*coca1*coca1+vyyjt*coca2*coca2+vzzjt*coca3*coca3
      !    1  - 2.0*(vxyjt*coca1*coca2+vxzjt*coca1*coca3+vyzjt*coca2*coca3)
      !       wxxj = wxxjt*coca1*coca1+wyyjt*coca2*coca2+wzzjt*coca3*coca3
      !    1  - 2.0*(wxyjt*coca1*coca2+wxzjt*coca1*coca3+wyzjt*coca2*coca3)

      uxxj = (sibb+sicc)*unji*coca1*coca1 - sibb*unji*cocb1*cocb1 &
            - sicc*unji*cocc1*cocc1 + 2.d0*(dft*coca1*cocb1 &
            + ddft*coca1*cocc1 - unji*sibc*cocb1*cocc1)+fuju*coca1*coca1
      vxxj = (sibb+sicc)*vnji*coca1*coca1 - sibb*vnji*cocb1*cocb1 &
            - sicc*vnji*cocc1*cocc1 + 2.d0*(dfn*coca1*cocb1 &
            + ddfn*coca1*cocc1 - vnji*sibc*cocb1*cocc1)+fvjv*coca1*coca1
      wxxj = (sibb+sicc)*wnji*coca1*coca1 - sibb*wnji*cocb1*cocb1 &
            - sicc*wnji*cocc1*cocc1 + 2.d0*(dfm*coca1*cocb1 &
            + ddfm*coca1*cocc1 - wnji*sibc*cocb1*cocc1)+fwjw*coca1*coca1

      !       print *, "xfix", sibb,sicc,unji,unbj,uncj,coca1,cocb1,cocc1
      !       print *, "xfix", uxxj, vxxj, wxxj
      !       stop

      !       if (abs(curv) > 5.0) then
   !  if (abs(curv) > 0.5d0/hx) then
   !     uxxj = 0.d0; vxxj = 0.d0; wxxj = 0.d0
   !  end if

      ucor = xcod(xg,x1,0.0d0,uxj,uxxj)
      vcor = xcod(xg,x1,0.0d0,vxj,vxxj)
      wcor = xcod(xg,x1,0.0d0,wxj,wxxj)

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine xfix here]
   subroutine xfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
         ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi,x,y,z,xi,u,v,w &
         ,ui,vi,wi,u1,v1,w1,cinfox &
         ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)

      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),vi(0:m,0:n,0:l),wi(0:m,0:n,0:l) &
            ,cinfox(np*4*n2),cinfo(np*n2) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,ft(n2),fn(n2),fm(n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex,kdm
      _REAL_ xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ddfm,ddfn,ddft
      _REAL_ dfm,dfn,dft,dsi,dsi2,fmv,fnv,ftv,hx1,hx2,sibb,sicc,sibc
      _REAL_ xg,yi,zi,fuju,fvjv,fwjw,ucor,vcor,wcor  
         
  !    ! common /radius/r0
      infodex=1
      hx2 = 2.0*hx
      hx1 = hx*hx

      yi = y(j)
      zi = z(k)
      if(ilr == 1) then
         xg = x(i+1)
      else
         xg = x(i-1)
      end if

      coca1 = cinfox(nx*np-np+4)
      coca2 = cinfox(nx*np-np+5)
      coca3 = cinfox(nx*np-np+6)
      cocb1 = cinfox(nx*np-np+7)
      cocb2 = cinfox(nx*np-np+8)
      cocb3 = cinfox(nx*np-np+9)
      cocc1 = cinfox(nx*np-np+10)
      cocc2 = cinfox(nx*np-np+11)
      cocc3 = cinfox(nx*np-np+12)
      sibb = cinfox(nx*np-np+13)
      sicc = cinfox(nx*np-np+14)
      sibc = cinfox(nx*np-np+15)
      curv = cinfox(nx*np-np+16)
      kdm = nint(cinfox(nx*np-np+22))
      !       if(master) then
      !         if (i==24.and.j==26.and.k==24) then
      !           print *, "x(24,26,24)",curv,sibb,sicc,sibc
      !         end if
      !       end if
      call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb &
            ,ft,fn,fm,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,curv,xi,yi,zi,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,infodex,nx,kdm,phi)
      
       
      !         if (master) then
      !         r2 = (xi**2+yi**2+zi**2)
      !         r1 = sqrt(r2)
      !         r4 = r2*r2
      !         cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !         print *, '111',unaj,unbj,uncj
      !         print *, '111',vnaj,vnbj,vncj
      !         print *, '111',wnaj,wnbj,wncj
      !         unx = -3.d0*cc*(yi*yi+zi*zi-3.d0*xi*xi)/r2/r4
      !         uny = -3.d0*cc*(-4.d0*xi*yi)/r2/r4
      !         unz = -3.d0*cc*(-4.d0*xi*zi)/r2/r4
      !         vnx = -3.d0*cc*(-4.d0*yi*xi)/r2/r4
      !         vny = -3.d0*cc*(xi*xi+zi*zi-3.d0*yi*yi)/r2/r4
      !         vnz = -3.d0*cc*(-4.d0*yi*zi)/r2/r4
      !         wnx = -3.d0*cc*(-4.d0*zi*xi)/r2/r4
      !         wny = -3.d0*cc*(-4.d0*zi*yi)/r2/r4
      !         wnz = -3.d0*cc*(xi*xi+yi*yi-3.d0*zi*zi)/r2/r4
      !         print *, '222',unx*coca1+uny*coca2+unz*coca3
      !    1      ,unx*cocb1+uny*cocb2+unz*cocb3
      !    1      ,unx*cocc1+uny*cocc2+unz*cocc3
      !         print *, '222',vnx*coca1+vny*coca2+vnz*coca3
      !    1      ,vnx*cocb1+vny*cocb2+vnz*cocb3
      !    1      ,vnx*cocc1+vny*cocc2+vnz*cocc3
      !         print *, '222',wnx*coca1+wny*coca2+wnz*coca3
      !    1      ,wnx*cocb1+wny*cocb2+wnz*cocb3
      !    1      ,wnx*cocc1+wny*cocc2+wnz*cocc3
      !         write(82,'(2f14.8)') ftv,-3.d0*cc*xi/(r4)
      !         write(82,'(2f14.8)') fnv,-3.d0*cc*yi/(r4)
      !         write(82,'(2f14.8)') fmv,-3.d0*cc*zi/(r4)
      !         print *,'fjump',xi,yi,zi,r2,-3.d0*cc*xi/r4,-3.d0*cc*yi/r4
      !    1      ,-3.d0*cc*zi/r4,ftv,fnv,fmv
      !         end if
      !         stop

      !        r2 = (xi**2+yi**2+zi**2)
      !        r1 = sqrt(r2)
      !        r4 = r2*r2
      !        cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !        unx = -3.d0*cc*(yi*yi+zi*zi-3.d0*xi*xi)/r2/r4
      !        uny = -3.d0*cc*(-4.d0*xi*yi)/r2/r4
      !        unz = -3.d0*cc*(-4.d0*xi*zi)/r2/r4
      !        vnx = -3.d0*cc*(-4.d0*yi*xi)/r2/r4
      !        vny = -3.d0*cc*(xi*xi+zi*zi-3.d0*yi*yi)/r2/r4
      !        vnz = -3.d0*cc*(-4.d0*yi*zi)/r2/r4
      !        wnx = -3.d0*cc*(-4.d0*zi*xi)/r2/r4
      !        wny = -3.d0*cc*(-4.d0*zi*yi)/r2/r4
      !        wnz = -3.d0*cc*(xi*xi+yi*yi-3.d0*zi*zi)/r2/r4
      !        unaj = unx*coca1+uny*coca2+unz*coca3
      !        unbj = unx*cocb1+uny*cocb2+unz*cocb3
      !        uncj = unx*cocc1+uny*cocc2+unz*cocc3
      !        vnaj = vnx*coca1+vny*coca2+vnz*coca3
      !        vnbj = vnx*cocb1+vny*cocb2+vnz*cocb3
      !        vncj = vnx*cocc1+vny*cocc2+vnz*cocc3
      !        wnaj = wnx*coca1+wny*coca2+wnz*coca3
      !        wnbj = wnx*cocb1+wny*cocb2+wnz*cocb3
      !        wncj = wnx*cocc1+wny*cocc2+wnz*cocc3
      !        ftv = -3.d0*cc*xi/(r4)
      !        fnv = -3.d0*cc*yi/(r4)
      !        fmv = -3.d0*cc*zi/(r4)
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. kdm) then
    ftv=ft(nn)
    fnv=fn(nn)
    fmv=fm(nn)
    endif
    endif       

      call xcor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,xg &
            ,xi,yi,zi,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor,ftv,dft,ddft &
            ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,fuju,fvjv,fwjw)

      !       call fjump(m,n,n2,np,i,j,index2,hx,hy,a,b,c,d,
      !    1         ft1,fn1,x,y,cinfo,coc,sic,curv,xi,yi,ftv,dft,ddft
      !    2         ,fnv,dfn,ddfn,ww1)

      !       ftv = -2.5*rnu*fti(t-dt)

      !       tb = t -dt

      !       call xcor(m,n,l,xa,xb,ya,yb,za,zb,tb,hx,hy,hz,phi,x,y,z,xg,
      !    1       xi,yi,zi,rnu,coc1,coc2,coc3,curv,fnv,dfn,ddfn,
      !    2       ftv,dft,ddft,fmv,dfm,ddfm,u,v,w,ucor1,vcor1,wcor1)

      !       call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb
      !    1    ,ft2,fn2,fm2,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,xi,yi,zi
      !    2    ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ww1,phi)

      !       ftv = -2.5*rnu*fti(t+dt)   !# should use the most updated (u,v)

      !       t1 = t + dt

      !       call xcor(m,n,l,xa,xb,ya,yb,za,zb,t1,hx,hy,hz,phi,x,y,z,xg
      !    1    ,xi,yi,zi,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      !    2    ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,ucor2,vcor2,wcor2)

      dsi = -float(ilr)
      dsi2 = float(info)

      !       if (i==46.and.j==47.and.k==29) print *,'inxfix',ui(i,j,k)
      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hx1
      ui(i,j,k)=ui(i,j,k) + dsi2*rnu*(ucor)/hx1
      !    1        + 0.0*dsi2*dsi*(1.5*u(i,j)*ucor/hx2-0.5*u1(i,j)*ucor1/hx2)

      !       if (i==18.and.j==31.and.k==20) print *,'xfix',ucor,ucor2,hx1
      !    1    ,rnu*ucor/hx1

      !       vi(i,j,k)=vi(i,j,k) + dsi2*rnu*0.5*(vcor+vcor2)/hx1
      vi(i,j,k)=vi(i,j,k) + dsi2*rnu*(vcor)/hx1
      !    1      +0.0*dsi2*dsi*(1.5*u(i,j)*vcor/hx2-0.5*u1(i,j)*vcor1/hx2)

      !       wi(i,j,k)=wi(i,j,k) + dsi2*rnu*0.5*(wcor+wcor2)/hx1
      wi(i,j,k)=wi(i,j,k) + dsi2*rnu*(wcor)/hx1
      !       p1(i,j) = p1(i,j) + dsi2*dsi*ucor/hx2
      !       if (i==46.and.j==47.and.k==29) print *,'inxfix',ui(i,j,k)
      !    1    ,ucor,ucor2

      return
   end subroutine

   !# ------------------------------------------------------------------

   double precision function xcod(xg,x1,uj,uxj,uxxj)

   implicit none
   _REAL_ xcod,x1,uj,uxj,uxxj,xg 
    
   if(solverorder .eq. 1) then
    
   xcod = uj + uxj*(xg - x1) 
   else
   xcod = uj + uxj*(xg - x1) + uxxj*(xg - x1)*(xg - x1)/2.0d0
   endif 
   return
end function
!FIXME:                nse.f, line 9128: bad indentation level at end of subroutine xfix




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ycor here]
   subroutine ycor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z &
         ,yg,x1,y1,z1,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor&
         ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w&
         ,fuju,fvjv,fwjw)

      implicit none
       _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l)

      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,x1,y1,z1,rnu
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,unbj,vnbj,wnbj,uncj,vncj,wncj
      _REAL_ ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ucor,vcor,wcor,unji
      _REAL_ uyj,uyyj,vnji,vyj,vyyj,wnji,wyj,wyyj,yg,fuju,fvjv,fwjw
      !       call fuj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fuju)
      !       call fvj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fvjv)
      !       call fwj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fwjw)

      unji = ftv
      vnji = fnv
      wnji = fmv

      !       uxj = unji*coca1
      !       vxj = vnji*coca1
      !       wxj = wnji*coca1
      uyj = unji*coca2
      vyj = vnji*coca2
      wyj = wnji*coca2
      !       uzj = unji*coca3
      !       vzj = vnji*coca3
      !       wzj = wnji*coca3

      !       uyyjt = -unji*curv    !# = (ft * tau *curv )/rnu
      !       vyyjt = -vnji*curv
      !       wyyjt = -wnji*curv

      !       uxyjt = dft
      !       vxyjt = dfn
      !       wxyjt = dfm

      !       uxzjt = dft
      !       vxzjt = dfn
      !       wxzjt = dfm

      !       uyzjt = dft
      !       vyzjt = dfn
      !       wyzjt = dfm

      !       call bilinear_m(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1
      !    1    ,x,y,z,phi,u,v,w,uxj,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj
      !    2    ,u1p,v1p,w1p)

      !       uxxjt = -uyyjt + ( - fuju)/rnu
      !       vxxjt = -vyyjt + ( - fvjv)/rnu
      !       wxxjt = -wyyjt + ( - fwjw)/rnu

      !       uzzjt = -uyyjt + ( - fuju)/rnu
      !       vzzjt = -vyyjt + ( - fvjv)/rnu
      !       wzzjt = -wyyjt + ( - fwjw)/rnu

      !       call gettan(coc1,coc2,coc3,coca1,coca2,coca3
      !    1    ,cocb1,cocb2,cocb3)

      !       coca = cocb1
      !       cocb = cocb2
      !       cocc = cocb3

      !       uyyj = uxxjt*coca*coca + uyyjt*cocb*cocb + uzzjt*cocc*cocc
      !    1  - 2.0*(uxyjt*coca*cocb + uxzjt*coca*cocc + uyzjt*cocb*cocc)
      !       vyyj = vxxjt*coca*coca + vyyjt*cocb*cocb + vzzjt*cocc*cocc
      !    1  - 2.0*(vxyjt*coca*cocb + vxzjt*coca*cocc + vyzjt*cocb*cocc)
      !       wyyj = wxxjt*coca*coca + wyyjt*cocb*cocb + wzzjt*cocc*cocc
      !    1  - 2.0*(wxyjt*coca*cocb + wxzjt*coca*cocc + wyzjt*cocb*cocc)

      uyyj = (sibb+sicc)*unji*coca2*coca2 - sibb*unji*cocb2*cocb2 &
            - sicc*unji*cocc2*cocc2 + 2.d0*(dft*coca2*cocb2 &
            + ddft*coca2*cocc2 - unji*sibc*cocb2*cocc2)+fuju*coca2*coca2
      vyyj = (sibb+sicc)*vnji*coca2*coca2 - sibb*vnji*cocb2*cocb2 &
            - sicc*vnji*cocc2*cocc2 + 2.d0*(dfn*coca2*cocb2 &
            + ddfn*coca2*cocc2 - vnji*sibc*cocb2*cocc2)+fvjv*coca2*coca2
      wyyj = (sibb+sicc)*wnji*coca2*coca2 - sibb*wnji*cocb2*cocb2 &
            - sicc*wnji*cocc2*cocc2 + 2.d0*(dfm*coca2*cocb2 &
            + ddfm*coca2*cocc2 - wnji*sibc*cocb2*cocc2)+fwjw*coca2*coca2

      !       if (abs(curv) > 5.0) then
  !   if (abs(curv) > 0.5d0/hx) then
  !      uyyj = 0.d0; vyyj = 0.d0; wyyj = 0.d0
  !   end if

      ucor = xcod(yg,y1,0.0d0,uyj,uyyj)
      vcor = xcod(yg,y1,0.0d0,vyj,vyyj)
      wcor = xcod(yg,y1,0.0d0,wyj,wyyj)

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine yfix here]
   subroutine yfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
         ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi,x,y,z,yj,u,v,w &
         ,ui,vi,wi,u1,v1,w1,cinfoy &
         ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)


      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),vi(0:m,0:n,0:l),wi(0:m,0:n,0:l) &
            ,cinfoy(np*4*n2),cinfo(np*n2) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,ft(n2),fn(n2),fm(n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex,kdm
      _REAL_ xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ddfm,ddfn,ddft
      _REAL_ dfm,dfn,dft,dsi,dsi2,fmv,fnv,ftv,hy1,hy2,sibb,sicc,sibc
      _REAL_ yg,yj,zi,fuju,fvjv,fwjw,ucor,vcor,wcor
      ! common /radius/r0

      hy2 = 2.0*hy
      hy1 = hy*hy
      infodex=2
      xi = x(i)
      zi = z(k)
      if(ilr == 1) then
         yg = y(j+1)
      else
         yg = y(j-1)
      end if

      coca1 = cinfoy(nx*np-np+4)
      coca2 = cinfoy(nx*np-np+5)
      coca3 = cinfoy(nx*np-np+6)
      cocb1 = cinfoy(nx*np-np+7)
      cocb2 = cinfoy(nx*np-np+8)
      cocb3 = cinfoy(nx*np-np+9)
      cocc1 = cinfoy(nx*np-np+10)
      cocc2 = cinfoy(nx*np-np+11)
      cocc3 = cinfoy(nx*np-np+12)
      sibb = cinfoy(nx*np-np+13)
      sicc = cinfoy(nx*np-np+14)
      sibc = cinfoy(nx*np-np+15)
      curv = cinfoy(nx*np-np+16)
      kdm=nint(cinfoy(nx*np-np+22))
      !       if(master) then
      !         if (i==24.and.j==26.and.k==24) then
      !           print *, "y(24,26,24)",curv,sibb,sicc,sibc
      !         end if
      !       end if
      call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb &
            ,ft,fn,fm,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,curv,xi,yj,zi &
            ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,infodex,nx,kdm,phi)
      
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. kdm) then
     ftv=ft(nn)
     fnv=fn(nn)
     fmv=fm(nn)
     endif
     endif
      !         r2 = (xi**2+yj**2+zi**2)
      !         r1 = sqrt(r2)
      !         r4 = r2*r2
      !         cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !         print *, '333/,unaj,unbj,uncj
      !         print *, '333',vnaj,vnbj,vncj
      !         print *, '333',wnaj,wnbj,wncj
      !         unx = -3.d0*cc*(yj*yj+zi*zi-3.d0*xi*xi)/r2/r4
      !         uny = -3.d0*cc*(-4.d0*xi*yj)/r2/r4
      !         unz = -3.d0*cc*(-4.d0*xi*zi)/r2/r4
      !         vnx = -3.d0*cc*(-4.d0*yj*xi)/r2/r4
      !         vny = -3.d0*cc*(xi*xi+zi*zi-3.d0*yj*yj)/r2/r4
      !         vnz = -3.d0*cc*(-4.d0*yj*zi)/r2/r4
      !         wnx = -3.d0*cc*(-4.d0*zi*xi)/r2/r4
      !         wny = -3.d0*cc*(-4.d0*zi*yj)/r2/r4
      !         wnz = -3.d0*cc*(xi*xi+yj*yj-3.d0*zi*zi)/r2/r4
      !         print *, '444',unx*coca1+uny*coca2+unz*coca3
      !    1      ,unx*cocb1+uny*cocb2+unz*cocb3
      !    1      ,unx*cocc1+uny*cocc2+unz*cocc3
      !         print *, '444',vnx*coca1+vny*coca2+vnz*coca3
      !    1      ,vnx*cocb1+vny*cocb2+vnz*cocb3
      !    1      ,vnx*cocc1+vny*cocc2+vnz*cocc3
      !         print *, '444',wnx*coca1+wny*coca2+wnz*coca3
      !    1      ,wnx*cocb1+wny*cocb2+wnz*cocb3
      !    1      ,wnx*cocc1+wny*cocc2+wnz*cocc3

      !        r2 = (xi**2+yj**2+zi**2)
      !        r1 = sqrt(r2)
      !        r4 = r2*r2
      !        cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !        unx = -3.d0*cc*(yj*yj+zi*zi-3.d0*xi*xi)/r2/r4
      !        uny = -3.d0*cc*(-4.d0*xi*yj)/r2/r4
      !        unz = -3.d0*cc*(-4.d0*xi*zi)/r2/r4
      !        vnx = -3.d0*cc*(-4.d0*yj*xi)/r2/r4
      !        vny = -3.d0*cc*(xi*xi+zi*zi-3.d0*yj*yj)/r2/r4
      !        vnz = -3.d0*cc*(-4.d0*yj*zi)/r2/r4
      !        wnx = -3.d0*cc*(-4.d0*zi*xi)/r2/r4
      !        wny = -3.d0*cc*(-4.d0*zi*yj)/r2/r4
      !        wnz = -3.d0*cc*(xi*xi+yj*yj-3.d0*zi*zi)/r2/r4
      !        unaj = unx*coca1+uny*coca2+unz*coca3
      !        unbj = unx*cocb1+uny*cocb2+unz*cocb3
      !        uncj = unx*cocc1+uny*cocc2+unz*cocc3
      !        vnaj = vnx*coca1+vny*coca2+vnz*coca3
      !        vnbj = vnx*cocb1+vny*cocb2+vnz*cocb3
      !        vncj = vnx*cocc1+vny*cocc2+vnz*cocc3
      !        wnaj = wnx*coca1+wny*coca2+wnz*coca3
      !        wnbj = wnx*cocb1+wny*cocb2+wnz*cocb3
      !        wncj = wnx*cocc1+wny*cocc2+wnz*cocc3
      !        ftv = -3.d0*cc*xi/(r4)
      !        fnv = -3.d0*cc*yj/(r4)
      !        fmv = -3.d0*cc*zi/(r4)

      call ycor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,yg &
            ,xi,yj,zi,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor,ftv,dft,ddft &
            ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,fuju,fvjv,fwjw)

      !       call fjump(m,n,n2,np,i,j,index2,hx,hy,a,b,c,d,
      !    1         ft1,fn1,x,y,cinfo,coc,sic,curv,xi,yj,ftv,dft,ddft
      !    2         ,fnv,dfn,ddfn,ww1)

      !       ftv = -2.5*rnu*fti(t-dt)

      !       tb = t -dt

      !       call xcor(m,n,l,xa,xb,ya,yb,za,zb,tb,hx,hy,hz,phi,x,y,z,xg,
      !    1       xi,yj,zi,rnu,coc1,coc2,coc3,curv,fnv,dfn,ddfn,
      !    2       ftv,dft,ddft,fmv,dfm,ddfm,u,v,w,ucor1,vcor1,wcor1)

      !       call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb
      !    1    ,ft2,fn2,fm2,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,xi,yj,zi
      !    2    ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ww1,phi)

      !       ftv = -2.5*rnu*fti(t+dt)   !# should use the most updated (u,v)

      !       t1 = t + dt

      !       call ycor(m,n,l,xa,xb,ya,yb,za,zb,t1,hx,hy,hz,phi,x,y,z,yg
      !    1    ,xi,yj,zi,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      !    2    ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,ucor2,vcor2,wcor2)

      dsi = -float(ilr)
      dsi2 = float(info)

      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hy1
      ui(i,j,k)=ui(i,j,k) + dsi2*rnu*(ucor)/hy1
      !    1        + 0.0*dsi2*dsi*(1.5*u(i,j)*ucor/hx2-0.5*u1(i,j)*ucor1/hx2)

      !       vi(i,j,k)=vi(i,j,k) + dsi2*rnu*0.5*(vcor+vcor2)/hy1
      vi(i,j,k)=vi(i,j,k) + dsi2*rnu*(vcor)/hy1
      !    1      +0.0*dsi2*dsi*(1.5*u(i,j)*vcor/hx2-0.5*u1(i,j)*vcor1/hx2)

      !       wi(i,j,k)=wi(i,j,k) + dsi2*rnu*0.5*(wcor+wcor2)/hy1
      wi(i,j,k)=wi(i,j,k) + dsi2*rnu*(wcor)/hy1
      !       p1(i,j) = p1(i,j) + dsi2*dsi*ucor/hx2

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zcor here]
   subroutine zcor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z &
         ,zg,x1,y1,z1,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
         ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor &
         ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,fuju,fvjv,fwjw)

      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l)
      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,x1,y1,z1,rnu
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,unbj,vnbj,wnbj,uncj,vncj,wncj
      _REAL_ ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ucor,vcor,wcor,unji
      _REAL_ uzj,uzzj,vnji,vzj,vzzj,wnji,wzj,wzzj,zg,fuju,fvjv,fwjw

      !       call fuj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fuju)
      !       call fvj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fvjv)
      !       call fwj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fwjw)

      unji = ftv
      vnji = fnv
      wnji = fmv

      !       uxj = unji*coca1
      !       vxj = vnji*coca1
      !       wxj = wnji*coca1
      !       uyj = unji*coca2
      !       vyj = vnji*coca2
      !       wyj = wnji*coca2
      uzj = unji*coca3
      vzj = vnji*coca3
      wzj = wnji*coca3

      !       uyyjt = -unji*curv    !# = (ft * tau *curv )/rnu
      !       vyyjt = -vnji*curv
      !       wyyjt = -wnji*curv

      !       uxyjt = dft
      !       vxyjt = dfn
      !       wxyjt = dfm

      !       uxzjt = dft
      !       vxzjt = dfn
      !       wxzjt = dfm

      !       uyzjt = dft
      !       vyzjt = dfn
      !       wyzjt = dfm

      !       call bilinear_m(m,n,l,xa,xb,ya,yb,za,zb,hx,hy,hz,x1,y1,z1
      !    1    ,x,y,z,phi,u,v,w,uxj,uyj,uzj,vxj,vyj,vzj,wxj,wyj,wzj
      !    2    ,u1p,v1p,w1p)

      !       uxxjt = -uyyjt + ( - fuju)/rnu
      !       vxxjt = -vyyjt + ( - fvjv)/rnu
      !       wxxjt = -wyyjt + ( - fwjw)/rnu

      !       uzzjt = -uyyjt + ( - fuju)/rnu
      !       vzzjt = -vyyjt + ( - fvjv)/rnu
      !       wzzjt = -wyyjt + ( - fwjw)/rnu

      !       call gettan(coc1,coc2,coc3,coca1,coca2,coca3
      !    1    ,cocb1,cocb2,cocb3)

      !       coca = cocc1
      !       cocb = cocc2
      !       cocc = cocc3

      !       uyyj = uxxjt*coca*coca + uyyjt*cocb*cocb + uzzjt*cocc*cocc
      !    1  - 2.0*(uxyjt*coca*cocb + uxzjt*coca*cocc + uyzjt*cocb*cocc)
      !       vyyj = vxxjt*coca*coca + vyyjt*cocb*cocb + vzzjt*cocc*cocc
      !    1  - 2.0*(vxyjt*coca*cocb + vxzjt*coca*cocc + vyzjt*cocb*cocc)
      !       wyyj = wxxjt*coca*coca + wyyjt*cocb*cocb + wzzjt*cocc*cocc
      !    1  - 2.0*(wxyjt*coca*cocb + wxzjt*coca*cocc + wyzjt*cocb*cocc)

      uzzj = (sibb+sicc)*unji*coca3*coca3 - sibb*unji*cocb3*cocb3 &
            - sicc*unji*cocc3*cocc3 + 2.d0*(dft*coca3*cocb3 &
            + ddft*coca3*cocc3 - unji*sibc*cocb3*cocc3)+fuju*coca3*coca3
      vzzj = (sibb+sicc)*vnji*coca3*coca3 - sibb*vnji*cocb3*cocb3 &
            - sicc*vnji*cocc3*cocc3 + 2.d0*(dfn*coca3*cocb3 &
            + ddfn*coca3*cocc3 - vnji*sibc*cocb3*cocc3)+fvjv*coca3*coca3
      wzzj = (sibb+sicc)*wnji*coca3*coca3 - sibb*wnji*cocb3*cocb3 &
            - sicc*wnji*cocc3*cocc3 + 2.d0*(dfm*coca3*cocb3 &
            + ddfm*coca3*cocc3 - wnji*sibc*cocb3*cocc3)+fwjw*coca3*coca3

      !       if (abs(curv) > 5.0) then
  !   if (abs(curv) > 0.5d0/hx) then
  !      uzzj = 0.d0; vzzj = 0.d0; wzzj = 0.d0
  !   end if

      ucor = xcod(zg,z1,0.0d0,uzj,uzzj)
      vcor = xcod(zg,z1,0.0d0,vzj,vzzj)
      wcor = xcod(zg,z1,0.0d0,wzj,wzzj)

      return
   end subroutine


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zfix here]
   subroutine zfix(m,n,l,nx,np,n2,i,j,k,info,ilr &
         ,xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,phi,x,y,z,zk,u,v,w &
         ,ui,vi,wi,u1,v1,w1,cinfoz &
         ,index2,cinfo,ft,fn,fm,fuju,fvjv,fwjw)


      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),vi(0:m,0:n,0:l),wi(0:m,0:n,0:l) &
            ,cinfoz(np*4*n2),cinfo(np*n2) &
            ,u1(0:m,0:n,0:l),v1(0:m,0:n,0:l),w1(0:m,0:n,0:l) &
            ,u(0:m,0:n,0:l),v(0:m,0:n,0:l),w(0:m,0:n,0:l) &
            ,ft(n2),fn(n2),fm(n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex,kdm
      _REAL_ xa,xb,ya,yb,za,zb,rnu,dt,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ddfm,ddfn,ddft
      _REAL_ dfm,dfn,dft,dsi,dsi2,fmv,fnv,ftv,hz1,hz2,sibb,sicc,sibc
      _REAL_ zg,zk,yi,fuju,fvjv,fwjw,ucor,vcor,wcor 
      ! common /radius/r0

      hz2 = 2.0*hz
      hz1 = hz*hz
     infodex=3
      xi = x(i)
      yi = y(j)
      if(ilr == 1) then
         zg = z(k+1)
      else
         zg = z(k-1)
      end if

      coca1 = cinfoz(nx*np-np+4)
      coca2 = cinfoz(nx*np-np+5)
      coca3 = cinfoz(nx*np-np+6)
      cocb1 = cinfoz(nx*np-np+7)
      cocb2 = cinfoz(nx*np-np+8)
      cocb3 = cinfoz(nx*np-np+9)
      cocc1 = cinfoz(nx*np-np+10)
      cocc2 = cinfoz(nx*np-np+11)
      cocc3 = cinfoz(nx*np-np+12)
      sibb = cinfoz(nx*np-np+13)
      sicc = cinfoz(nx*np-np+14)
      sibc = cinfoz(nx*np-np+15)
      curv = cinfoz(nx*np-np+16)
      kdm = nint(cinfoz(nx*np-np+22))
      !       if(master) then
      !         if (i==24.and.j==26.and.k==24) then
      !           print *, "z(24,26,24)",curv,sibb,sicc,sibc
      !         end if
      !       end if 
      call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb &
            ,ft,fn,fm,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,curv,xi,yi,zk &
            ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,infodex,nx,kdm,phi)
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. kdm) then
     ftv=ft(nn)
     fnv=fn(nn)
     fmv=fm(nn)
     endif
     endif
      !         r2 = (xi**2+yi**2+zk**2)
      !         r1 = sqrt(r2)
      !         r4 = r2*r2
      !         cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !         print *, '555',unaj,unbj,uncj
      !         print *, '555',vnaj,vnbj,vncj
      !         print *, '555',wnaj,wnbj,wncj
      !         unx = -3.d0*cc*(yi*yi+zk*zk-3.d0*xi*xi)/r2/r4
      !         uny = -3.d0*cc*(-4.d0*xi*yi)/r2/r4
      !         unz = -3.d0*cc*(-4.d0*xi*zk)/r2/r4
      !         vnx = -3.d0*cc*(-4.d0*yi*xi)/r2/r4
      !         vny = -3.d0*cc*(xi*xi+zk*zk-3.d0*yi*yi)/r2/r4
      !         vnz = -3.d0*cc*(-4.d0*yi*zk)/r2/r4
      !         wnx = -3.d0*cc*(-4.d0*zk*xi)/r2/r4
      !         wny = -3.d0*cc*(-4.d0*zk*yi)/r2/r4
      !         wnz = -3.d0*cc*(xi*xi+yi*yi-3.d0*zk*zk)/r2/r4
      !         print *, '666',unx*coca1+uny*coca2+unz*coca3
      !    1      ,unx*cocb1+uny*cocb2+unz*cocb3
      !    1      ,unx*cocc1+uny*cocc2+unz*cocc3
      !         print *, '666',vnx*coca1+vny*coca2+vnz*coca3
      !    1      ,vnx*cocb1+vny*cocb2+vnz*cocb3
      !    1      ,vnx*cocc1+vny*cocc2+vnz*cocc3
      !         print *, '666',wnx*coca1+wny*coca2+wnz*coca3
      !    1      ,wnx*cocb1+wny*cocb2+wnz*cocb3
      !    1      ,wnx*cocc1+wny*cocc2+wnz*cocc3

      !        r2 = (xi**2+yi**2+zk**2)
      !        r1 = sqrt(r2)
      !        r4 = r2*r2
      !        cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !        unx = -3.d0*cc*(yi*yi+zk*zk-3.d0*xi*xi)/r2/r4
      !        uny = -3.d0*cc*(-4.d0*xi*yi)/r2/r4
      !        unz = -3.d0*cc*(-4.d0*xi*zk)/r2/r4
      !        vnx = -3.d0*cc*(-4.d0*yi*xi)/r2/r4
      !        vny = -3.d0*cc*(xi*xi+zk*zk-3.d0*yi*yi)/r2/r4
      !        vnz = -3.d0*cc*(-4.d0*yi*zk)/r2/r4
      !        wnx = -3.d0*cc*(-4.d0*zk*xi)/r2/r4
      !        wny = -3.d0*cc*(-4.d0*zk*yi)/r2/r4
      !        wnz = -3.d0*cc*(xi*xi+yi*yi-3.d0*zk*zk)/r2/r4
      !        unaj = unx*coca1+uny*coca2+unz*coca3
      !        unbj = unx*cocb1+uny*cocb2+unz*cocb3
      !        uncj = unx*cocc1+uny*cocc2+unz*cocc3
      !        vnaj = vnx*coca1+vny*coca2+vnz*coca3
      !        vnbj = vnx*cocb1+vny*cocb2+vnz*cocb3
      !        vncj = vnx*cocc1+vny*cocc2+vnz*cocc3
      !        wnaj = wnx*coca1+wny*coca2+wnz*coca3
      !        wnbj = wnx*cocb1+wny*cocb2+wnz*cocb3
      !        wncj = wnx*cocc1+wny*cocc2+wnz*cocc3
      !        ftv = -3.d0*cc*xi/(r4)
      !        fnv = -3.d0*cc*yi/(r4)
      !        fmv = -3.d0*cc*zk/(r4)

      call zcor(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,zg &
            ,xi,yi,zk,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,vcor,wcor,ftv,dft,ddft &
            ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,fuju,fvjv,fwjw)
      !       if (i==31.and.j==33.and.k==13) then
      !         r2 = (xi**2+yi**2+zk**2)
      !         r4 = r2*r2
      !         r0 = 1.2d0
      !         cc = 4.0*r0*rnu-sqrt(16.0*rnu**2*r0**2-2.0*(r0-r0**3)+1.d-16)
      !         print *,'fjump',xi,yi,zk,r2,-3.d0*cc*xi/r4,-3.d0*cc*yi/r4
      !    1      ,-3.d0*cc*zk/r4,ftv,fnv,fmv
      !         stop
      !       end if

      !       call fjump(m,n,n2,np,i,j,index2,hx,hy,a,b,c,d,
      !    1         ft1,fn1,x,y,cinfo,coc,sic,curv,xi,yi,ftv,dft,ddft
      !    2         ,fnv,dfn,ddfn,ww1)

      !       ftv = -2.5*rnu*fti(t-dt)

      !       tb = t -dt

      !       call xcor(m,n,l,xa,xb,ya,yb,za,zb,tb,hx,hy,hz,phi,x,y,z,xg,
      !    1       xi,yi,zk,rnu,coc1,coc2,coc3,curv,fnv,dfn,ddfn,
      !    2       ftv,dft,ddft,fmv,dfm,ddfm,u,v,w,ucor1,vcor1,wcor1)

      !       call fjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz,xa,xb,ya,yb,za,zb
      !    1    ,ft2,fn2,fm2,x,y,z,cinfo,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,xi,yi,zk
      !    2    ,ftv,dft,ddft,fnv,dfn,ddfn,fmv,dfm,ddfm,ww1,phi)

      !       ftv = -2.5*rnu*fti(t+dt)   !# should use the most updated (u,v)

      !       t1 = t + dt

      !       call zcor(m,n,l,xa,xb,ya,yb,za,zb,t1,hx,hy,hz,phi,x,y,z,zg
      !    1    ,xi,yi,zk,rnu,coca1,coca2,coca3,cocb1,cocb2,cocb3
      !    1    ,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      !    2    ,fnv,dfn,ddfn,fmv,dfm,ddfm,u,v,w,ucor2,vcor2,wcor2)

      !       if (i==31.and.j==33.and.k==13) print *,'zfix',ucor,hz1
      !    1    ,rnu*ucor/hz1

      dsi = -float(ilr)
      dsi2 = float(info)

      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hz1
      ui(i,j,k)=ui(i,j,k) + dsi2*rnu*(ucor)/hz1
      !    1        + 0.0*dsi2*dsi*(1.5*u(i,j)*ucor/hx2-0.5*u1(i,j)*ucor1/hx2)

      !       vi(i,j,k)=vi(i,j,k) + dsi2*rnu*0.5*(vcor+vcor2)/hz1
      vi(i,j,k)=vi(i,j,k) + dsi2*rnu*(vcor)/hz1
      !    1      +0.0*dsi2*dsi*(1.5*u(i,j)*vcor/hx2-0.5*u1(i,j)*vcor1/hx2)

      !       wi(i,j,k)=wi(i,j,k) + dsi2*rnu*0.5*(wcor+wcor2)/hz1
      wi(i,j,k)=wi(i,j,k) + dsi2*rnu*(wcor)/hz1
      !       p1(i,j) = p1(i,j) + dsi2*dsi*ucor/hx2

      return
   end subroutine
  
    INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  NUMROC computes the NUMber of Rows Or Columns of a distributed
!  matrix owned by the process indicated by IPROC.
!
!  Arguments
!  =========
!
!  N         (global input) INTEGER
!            The number of rows/columns in distributed matrix.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local input) INTEGER
!            The coordinate of the process whose local array row or
!            column is to be determined.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row or column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER              EXTRABLKS, MYDIST, NBLOCKS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
!     Figure PROC's distance from source process
!
      MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
!
!     Figure the total number of whole NB blocks N is split up into
!
      NBLOCKS = N / NB
!
!     Figure the minimum number of rows/cols a process can have
!
      NUMROC = (NBLOCKS/NPROCS) * NB
!
!     See if there are any extra blocks
!
      EXTRABLKS = MOD( NBLOCKS, NPROCS )
!
!     If I have an extra block
!
      IF( MYDIST.LT.EXTRABLKS ) THEN
          NUMROC = NUMROC + NB
!
!         If I have last block, it may be a partial block
!
      ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
          NUMROC = NUMROC + MOD( N, NB )
      END IF
!
      RETURN
!
!     End of NUMROC
!
END FUNCTION
       subroutine matvecp(m,n,l,np,n2,indexx,indexy,indexz &
          ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p,phi,cinfox,cinfoy,cinfoz&
          ,index2,cinfo,fn,fvec,bf,ubc2,bv)

      implicit none 

  !    ! common /radius/r0
  !    ! common /aret/aret
      !*************************************************************************
      !   This routine
      ! -----------------------------------------------------------------------
      ! x,y,z:coordinate      
      ! u,v,w: velocity
      ! p1,p0,p: pressure
      ! phi: level set
      ! index2: projection point
      ! cinfo: coordinate, normal,tangential of irregular points
      ! ww:working array
      !cinfox,cinfoy,cinfoz:coordinate,normal,tangential of irregular points
      !indexx,indexy,indexz:irregular points
      !ft,fn,fm:jump conditions um,un,ut
      !fec,ubc2:matrix vector multiplication result
      !bf:right-hand-side
      !aret:volume of each domain
      !ww1,ww2,vl,uw: SVD
      !kzone: domain of grid points 





      _REAL_ x(0:m),y(0:n),z(0:l)&
            ,p(0:m,0:n,0:l),bv(0:m,0:n,0:l)&
            ,phi(0:m,0:n,0:l)
     integer index2(0:m,0:n,0:l) 
     _REAL_  cinfo(np*n2) &
            ,cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2) 
     integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
     _REAL_ fn(n2),fvec(n2),ubc2(n2),bf(n2)
      integer i,j,k,l,m,n,np,n2
      _REAL_ xa,xb,ya,yb,za,zb,hx,hy,hz,t,rnu,bci1,bci2,bci3
      _REAL_ coca1,coca2,cc,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ curv,h,pjmp,e
      integer n22,nn
      _REAL_ usout,utout,uvout,uout 

      call pbc(m,n,l,np,n2,indexx,indexy,indexz &
            ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p &
            ,phi,cinfox,cinfoy,cinfoz,cinfo,index2,fn,bv)
         e=0.0d0

      do k=1,l-1
         do j=1,n-1
            do i=1,m-1
               nn = index2(i,j,k)
               n22 = abs(nn)
               ! was only outside, now inside --CQ
               if(nn < 0 ) then
                if(ifupdate .eqv. .true.) then
            call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                 cinfo(n22*np-np+22)=pjmp
                 else   
                 pjmp=cinfo(n22*np-np+22)
                 endif

         call interpp(m,n,l,np,n2,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,t,phi,fn,x,y,z &
         ,p,uvout,utout,usout,uout,pjmp,cinfo)

                  fvec(n22)=uout-ubc2(n22)
                  fvec(n22) = fvec(n22) + bf(n22)
       !           if(e<fvec(n22)) then
       !              e=fvec(n22)
       !            endif
   !        write(11,*) fn(n22),fvec(n22),p(i,j,k),ubc2(n22),bf(n22)
               end if  ! (nn < 0)
            end do  !  i=1,m-1
         end do  !  j=1,n-1
      end do  !  k=1,l-1
   !   write(9,*) e
      return
   end subroutine matvecp

      subroutine pbc(m,n,l,np,n2,indexx,indexy,indexz &
            ,xa,xb,ya,yb,za,zb,hx,hy,hz,t,x,y,z,rnu,p &
            ,phi,cinfox,cinfoy,cinfoz,cinfo,index2,fn,bv)
      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l) &
            ,ui(0:m,0:n,0:l),bv(0:m,0:n,0:l),p(0:m,0:n,0:l) &
            ,phi(0:m,0:n,0:l)
      integer index2(0:m,0:n,0:l)
      _REAL_  cinfox(np*4*n2),cinfoy(np*4*n2),cinfoz(np*4*n2),cinfo(np*n2) 
       integer indexx(0:m,0:n,0:l),indexy(0:m,0:n,0:l),indexz(0:m,0:n,0:l) 
       _REAL_ fn(n2)
      integer i,j,k,m,n,l,np,n2
      _REAL_ h,hx1,hx2,hy1,hy2,hz1,hz2,hx,hy,hz,rnu,xa,xb,zk
      _REAL_ ya,yb,za,zb,t,f3,fc,fl,pjmp
      integer ilr,info,nx,ny,nz,solvopt,optx,opty,optz
      _REAL_  x1,x2,x3,xi,yj
      _REAL_  epsx,epsy,epsz,iv(1:(m-1)*(n-1)*(l-1)),xso((m-1)*(n-1)*(l-1)+2*(m-1)*(n-1)),dummy,accept 
      ! common /radius/r0
      !-------------------------------initialize---------------------
      h = min(hx,hy,hz)
      hx1 = hx*hx
      hx2 = 2.0d0*hx
      hy1 = hy*hy
      hy2 = 2.0d0*hy
      hz1 = hz*hz
      hz2 = 2.0d0*hz
      solvopt=1
      if(ifzero .eqv. .true.) then
      p=0.0d0
      ifzero= .false.
      return
      endif
      do k=1,l-1
         do j=1,n-1
           do i=1,m-1
           ui(i,j,k)=bv(i,j,k)*h*h 
           enddo
        enddo
      enddo
!
      do  k=1,l-1
         do  j=1,n-1
            do  i=1,m-1
              ! x1 = x(i) - hx
               x2 = x(i)
               x3 = x(i) + hx
               fl = phi(i-1,j,k)
               fc = phi(i,j,k)
               f3 = phi(i+1,j,k)
               
               !# -------------------------------- !# Irregular grid?
               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexx(i,j,k)
                  xi = cinfox(nx*np-np+1)
                  !find correction terms
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfox(nx*np-np+21)
                  endif
                  call xfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,xi,ui,cinfox,cinfo,index2,fn,pjmp)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i-2,j,k) > 0.0) then
                     nx = indexx(i-1,j,k) + 1
                  else
                     nx = indexx(i-1,j,k)
                  end if
                  xi = cinfox(nx*np-np+1)
                  ilr = -1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfox(nx*np-np+21)
                  endif
                  call xfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,xi,ui,cinfox,cinfo,index2,fn,pjmp)
               end if                 !# End of (x_{i-1},x_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i-1,j,k) > 0.0) then
                     nx = indexx(i,j,k) + 1
                  else
                     nx = indexx(i,j,k)
                  end if
                  xi = cinfox(nx*np-np+1)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfox(nx*np-np+21)
                  endif
                  call xfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,xi,ui,cinfox,cinfo,index2,fn,pjmp)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexx(i+1,j,k)
                  xi = cinfox(nx*np-np+1)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfox(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfox(nx*np-np+21)
                  endif
                  call xfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,xi,ui,cinfox,cinfo,index2,fn,pjmp)
               end if                   !# End of (x_i,x_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
               !       if (i==34.and.j==30.and.k==13) print *,'uvw',ui(i,j,k)
               !    1    ,vi(i,j,k),wi(i,j,k)


               x1 = y(j) - hy
               x2 = y(j)
               x3 = y(j) + hy
               fl = phi(i,j-1,k)
               fc = phi(i,j,k)
               f3 = phi(i,j+1,k)
               !# ------------     !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexy(i,j,k)
                  yj = cinfoy(nx*np-np+2)
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoy(nx*np-np+21)
                  endif

                  call yfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,yj,ui,cinfoy,cinfo,index2,fn,pjmp)

               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i,j-2,k) > 0.0) then
                     nx = indexy(i,j-1,k) + 1
                  else
                     nx = indexy(i,j-1,k)
                  end if
                  yj = cinfoy(nx*np-np+2)
                  ilr = -1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoy(nx*np-np+21)
                  endif
                  call yfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,yj,ui,cinfoy,cinfo,index2,fn,pjmp)
               end if                 !# End of (x_{i-1},x_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i,j-1,k) > 0.0) then
                     nx = indexy(i,j,k) + 1
                  else
                     nx = indexy(i,j,k)
                  end if
                  yj = cinfoy(nx*np-np+2)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoy(nx*np-np+21)
                  endif
                  call yfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,yj,ui,cinfoy,cinfo,index2,fn,pjmp)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexy(i,j+1,k)
                  yj = cinfoy(nx*np-np+2)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoy(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoy(nx*np-np+21)
                  endif
                  call yfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,yj,ui,cinfoy,cinfo,index2,fn,pjmp)
               end if                   !# End of (x_i,x_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)


               x1 = z(k) - hz
               x2 = z(k)
               x3 = z(k) + hz
               fl = phi(i,j,k-1)
               fc = phi(i,j,k)
               f3 = phi(i,j,k+1)
               !# ------------     !# Irregular grid?

               if( fl > 0.0 .and. fc <= 0.0) then
                  !       if( fl .gt. 0.0 .and. fc .le. 0.0 .and. f3 .le. 0.0 ) then
                  info = 1
                  ilr = -1         !# dlr = 1, right, dlr =-1, left
                  nx = indexz(i,j,k)
                  zk = cinfoz(nx*np-np+3)
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoz(nx*np-np+21)
                  endif

                  call zfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,zk,ui,cinfoz,cinfo,index2,fn,pjmp)

               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( fl <= 0.0 .and. fc > 0.0) then
                  !       if( fl .le. 0.0 .and. fc .gt. 0.0 .and. f3 .gt. 0.0 ) then
                  if(phi(i,j,k-2) > 0.0) then
                     nx = indexz(i,j,k-1) + 1
                  else
                     nx = indexz(i,j,k-1)
                  end if
                  zk = cinfoz(nx*np-np+3)
                  ilr = -1
                  info = -1
                 
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoz(nx*np-np+21)
                  endif
                  call zfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,zk,ui,cinfoz,cinfo,index2,fn,pjmp)
               end if                 !# End of (z_{i-1},z_i) -------------------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 > 0.0 .and. fc <= 0.0) then
                  !       if( f3 .gt. 0.0 .and. fc .le. 0.0 .and. fl .le. 0.0 ) then
                  if(phi(i,j,k-1) > 0.0) then
                     nx = indexz(i,j,k) + 1
                  else
                     nx = indexz(i,j,k)
                  end if
                  zk = cinfoz(nx*np-np+3)
                  ilr = 1
                  info = 1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoz(nx*np-np+21)
                  endif
                  call zfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,zk,ui,cinfoz,cinfo,index2,fn,pjmp)
               end if
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)

               if( f3 <= 0.0 .and. fc > 0.0) then
                  !       if( f3 .le. 0.0 .and. fc .gt. 0.0 .and. fl .gt. 0.0 ) then
                  nx = indexz(i,j,k+1)
                  zk = cinfoz(nx*np-np+3)
                  ilr = 1
                  info = -1
                  if(ifupdate .eqv. .true.) then   
                  call fjmps(pjmp,phi,bv(1:m-1,1:n-1,1:l-1),i,j,k,m-1,n-1,l-1)
                  cinfoz(nx*np-np+21)=pjmp
                  else           
                  pjmp=cinfoz(nx*np-np+21)
                  endif
                  call zfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,zk,ui,cinfoz,cinfo,index2,fn,pjmp)
               end if                   !# End of (z_i,z_{i+1}) ----------
               !       if (i==16.and.j==14.and.k==9) print *,'fix',ui(i,j,k)
             enddo
           enddo
         enddo
      !       write(*,*) 'bbb',ui(1,1,1)

         if(ifupdate .eqv. .true.) then
         solvernorm=0.0d0
         do k=1,l-1
            do j=1,n-1
               do i=1,m-1
               if(index2(i,j,k) .lt. 0) then
               solvernorm=solvernorm+abs(ui(i,j,k))
               endif
               enddo
             enddo
         enddo
         solvernorm=solvernorm/n2*h
        if(master) then
          write(6,*) "norm of pressure solver",solvernorm
          flush(6)
         endif
         endif

 
         if (nse_debug .le. 8 .and. nse_debug .ne.0) then
         do j = 1, n-1; do i = 1, m-1
               ui(i,j,1 ) = ui(i,j,1 ) - pe(rnu,t,x(i),y(j),z(0))
               p(i,j,0)=0 
               ui(i,j,l-1) = ui(i,j,l-1 )- pe(rnu,t,x(i),y(j),z(l))
               p(i,j,l)=0
         end do; end do
           
         ! j=0 and ym+1 faces
           
         do k = 1, l-1; do i = 1, m-1
               ui(i,1,k) = ui(i,1,k) - pe(rnu,t,x(i),y(0),z(k))
               p(i,0,k)=0
               ui(i,n-1,k) = ui(i,n-1,k ) - pe(rnu,t,x(i),y(n),z(k))
               p(i,n,k)=0
         end do; end do
      
         ! i=0 and i=xm+1 faces
      
         do k = 1, l-1; do j = 1, n-1
               ui(1,j,k) = ui(1,j,k ) - pe(rnu,t,x(0),y(j),z(k))
               p(0,j,k)=0
               ui(m-1,j,k) = ui(m-1,j,k ) - pe(rnu,t,x(m),y(j),z(k))
               p(m,j,k)=0 
         end do; end do
        endif

      accept=1e-8
      dummy=0.0d0
      solvopt=1
      epsx=1.d0
      epsy=1.d0
      epsz=1.d0
     
        optx=0
        opty=0
        optz=0
            call init_param(l-1,m-1,n-1,(l-1)*(m-1),(l-1)*(m-1)*(n-1),10000,dummy,accept,&
                 0.0d0,1.0d0,h,1.9d0)

            call allocate_array(solvopt) 
            xso(:) = 0.d0
            call init_array(solvopt,epsx,epsy,epsz,-ui(1:m-1,1:n-1,1:l-1),iv,xso,optx,opty,optz)
            call pb_iccg(p(1:m-1,1:n-1,1:l-1),xso)  
          call deallocate_array(solvopt)
  !   do k=1,l-1
  !      do j=1,n-1
  !         do i=1,m-1
  !            write(123,*)p(i,j,k),pe(rnu,t,x(i),y(j),z(k)),(p(i,j,k)-pe(rnu,t,x(i),y(j),z(k)))/p(i,j,k)
  !            write(456,*) ui(i,j,k)
  !                      write(123789123,*) i,j,k,ui(i,j,k),ui(i,j,k)-(p(i-1,j,k)+p(i+1,j,k)+p(i,j+1,k)&
  !                     +p(i,j-1,k)+p(i,j,k+1)+p(i,j,k-1)-6*p(i,j,k)),&
  !                  1-(p(i-1,j,k)+p(i+1,j,k)+p(i,j+1,k) +p(i,j-1,k)+p(i,j,k+1)+p(i,j,k-1)-6*p(i,j,k))/(ui(i,j,k))

  !         end do
  !      end do
  !   end do
  !    stop    
      return
   end subroutine pbc

   subroutine pfjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,fn,x,y,z,uji,ujd,ujdd,infodex,nx)
     implicit none
     integer nl,ne
      parameter( nl = 16, ne=6)
     

      !# ********************************************************************
      !#   This routine evaluate the jump information [u], [u]', [u]'', at a
      !#   point (x1,y1) on the interface using an arc-length interpolation.
      !#
      !#   Input information:
      !#
      !#   Onput information:
      !#     uji,ujd,ujdd, are the jumps of [u], [u]', [u]''.
      !# -----------------------------------------------------------------------

      _REAL_ x(0:m),y(0:n),z(0:l),w3(nl),w3b(nl),w3c(nl) 
      integer index2(0:m,0:n,0:l)
      _REAL_  ss(nl),fn(n2) 
      integer ix(nl),iy(nl),iz(nl)
      _REAL_ hx,hy,hz,xa,xb,ya,yb,za,zb  
      integer m,n,l,n2,np,i,j,k,nn,nn1,nx
      _REAL_ uji,ujd,ujdd
      _REAL_ h,ts
      integer i1,j1,k1,i2,infodex,nsub

   !   ! common /radius/r0
      !       write(54,'(a,3i6)') 'fjump',i,j,k
      nn=nx+(infodex-1)*2*n2
      nsub=numfjump(nn)
      do i1=1,nsub
       ix(i1)= nint(fjumpinfo((nn-1)*npfjump+1,i1)) 
       iy(i1)= nint(fjumpinfo((nn-1)*npfjump+2,i1))
       iz(i1)= nint(fjumpinfo((nn-1)*npfjump+3,i1)) 
       ss(i1)= fjumpinfo((nn-1)*npfjump+4,i1) 
       w3(i1) =fjumpinfo((nn-1)*npfjump+5,i1) 
       w3b(i1)=fjumpinfo((nn-1)*npfjump+6,i1)
       w3c(i1)=fjumpinfo((nn-1)*npfjump+7,i1)
      end do
      !# ------------------ Summarize to get u_j' and u_j'' ------------------

      uji = 0.0d0
      ujd = 0.0d0
      ujdd = 0.0d0
      do  i2=1, nsub     !# Starting from the center
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)
         nn1 = abs(index2(i1,j1,k1))

         uji = uji + w3(i2)*fn(nn1)*ts
         ujd = ujd + w3b(i2)*fn(nn1)*ts
         ujdd = ujdd + w3c(i2)*fn(nn1)*ts
      enddo
 
      return
   end subroutine pfjump

   subroutine xcorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,xg &
            ,xi,yi,zi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)


      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,xi,yi,zi
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,ftv,dft,ddft,ucor,unji,uxj,uxxj,pjmp
 
      !       call fuj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fuju)
      !       call fvj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fvjv)
      !       call fwj(m,n,l,t,rnu,phi,x,y,z,x1,y1,z1,coc1,coc2,coc3
      !    1    ,curv,fwjw)
      ! above outputs fuju fvjv fwjw all zero --CQ

      unji = ftv

      uxj = unji*coca1

      uxxj = (sibb+sicc)*unji*coca1*coca1 - sibb*unji*cocb1*cocb1 &
            - sicc*unji*cocc1*cocc1 + 2.d0*(dft*coca1*cocb1 &
            + ddft*coca1*cocc1 - unji*sibc*cocb1*cocc1)+pjmp*coca1*coca1


      ucor = xcod(xg,xi,0.0d0,uxj,uxxj)

      return
   end subroutine xcorp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine xfix here]
      subroutine xfixp(m,n,l,nx,np,n2,i,j,k,info,ilr &
                        ,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi &
                        ,x,y,z,xi,ui,cinfox,cinfo,index2,fn,pjmp)

      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),cinfox(np*4*n2),fn(n2),cinfo(np*n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      _REAL_ dsi,dsi2,hx1,hx2,sibb,sicc,sibc,xg,yi,zi,ucor,pjmp  
         
  !    ! common /radius/r0

      hx2 = 2.0*hx
      hx1 = hx*hx
      infodex=1
      yi = y(j)
      zi = z(k)
      if(ilr == 1) then
         xg = x(i+1)
      else
         xg = x(i-1)
      end if

      coca1 = cinfox(nx*np-np+4)
      coca2 = cinfox(nx*np-np+5)
      coca3 = cinfox(nx*np-np+6)
      cocb1 = cinfox(nx*np-np+7)
      cocb2 = cinfox(nx*np-np+8)
      cocb3 = cinfox(nx*np-np+9)
      cocc1 = cinfox(nx*np-np+10)
      cocc2 = cinfox(nx*np-np+11)
      cocc3 = cinfox(nx*np-np+12)
      sibb = cinfox(nx*np-np+13)
      sicc = cinfox(nx*np-np+14)
      sibc = cinfox(nx*np-np+15)
      curv = cinfox(nx*np-np+16)
      
      call pfjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,fn,x,y,z,ftv,dft,ddft,infodex,nx)
       
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. nint(cinfox(nx*np-np+22))) then
     ftv=fn(nn)
     endif  
     endif    

      call xcorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,xg &
            ,xi,yi,zi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)

      dsi = -float(ilr)
      dsi2 = float(info)

      !       if (i==46.and.j==47.and.k==29) print *,'inxfix',ui(i,j,k)
      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hx1
      ui(i,j,k)=ui(i,j,k) + dsi2*ucor
      !    1        + 0.0*dsi2*dsi*(1.5*u(i,j)*ucor/hx2-0.5*u1(i,j)*ucor1/hx2)

      return
   end subroutine xfixp



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ycor here]
     subroutine ycorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,yg &
            ,xi,yi,zi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)

      implicit none
       _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 

      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,xi,yi,zi
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,ftv,dft,ddft,ucor,unji,uyj,uyyj,yg,pjmp

      unji = ftv
      uyj = unji*coca2

      uyyj = (sibb+sicc)*unji*coca2*coca2 - sibb*unji*cocb2*cocb2 &
            - sicc*unji*cocc2*cocc2 + 2.d0*(dft*coca2*cocb2 &
            + ddft*coca2*cocc2 - unji*sibc*cocb2*cocc2)+pjmp*coca2*coca2

      ucor = xcod(yg,yi,0.0d0,uyj,uyyj)

      return
   end subroutine ycorp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine yfix here]
       subroutine yfixp(m,n,l,nx,np,n2,i,j,k,info,ilr,xa,xb,ya,yb,za,zb&
           ,t,hx,hy,hz,phi,x,y,z,yj,ui,cinfoy,cinfo,index2,fn,pjmp)


      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),cinfoy(np*4*n2),fn(n2),cinfo(np*n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      _REAL_ dsi,dsi2,hy1,hy2,sibb,sicc,sibc
      _REAL_ yg,yj,zi,pjmp,ucor
      ! common /radius/r0

      hy2 = 2.0*hy
      hy1 = hy*hy
      infodex=2
      xi = x(i)
      zi = z(k)
      if(ilr == 1) then
         yg = y(j+1)
      else
         yg = y(j-1)
      end if

      coca1 = cinfoy(nx*np-np+4)
      coca2 = cinfoy(nx*np-np+5)
      coca3 = cinfoy(nx*np-np+6)
      cocb1 = cinfoy(nx*np-np+7)
      cocb2 = cinfoy(nx*np-np+8)
      cocb3 = cinfoy(nx*np-np+9)
      cocc1 = cinfoy(nx*np-np+10)
      cocc2 = cinfoy(nx*np-np+11)
      cocc3 = cinfoy(nx*np-np+12)
      sibb = cinfoy(nx*np-np+13)
      sicc = cinfoy(nx*np-np+14)
      sibc = cinfoy(nx*np-np+15)
      curv = cinfoy(nx*np-np+16)
      call pfjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,fn,x,y,z,ftv,dft,ddft,infodex,nx)
      
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. nint(cinfoy(nx*np-np+22))) then
      ftv=fn(nn)
      endif
     endif

      call ycorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,yg &
            ,xi,yj,zi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)


      dsi = -float(ilr)
      dsi2 = float(info)

      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hy1
      ui(i,j,k)=ui(i,j,k) + dsi2*ucor

      return
   end subroutine yfixp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zcor here]
     subroutine zcorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,zg &
            ,xi,yi,zi,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)

      implicit none
      _REAL_ x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer m,n,l
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xg,xi,yi,zi
      _REAL_ coca1,coca2,coca3,cocb1,cocb2,cocb3,cocc1,cocc2,cocc3
      _REAL_ sibb,sicc,sibc,curv,ftv,dft,ddft,ucor,unji,uzj,uzzj,zg,pjmp


      unji = ftv

      uzj = unji*coca3

      uzzj = (sibb+sicc)*unji*coca3*coca3 - sibb*unji*cocb3*cocb3 &
            - sicc*unji*cocc3*cocc3 + 2.d0*(dft*coca3*cocb3 &
            + ddft*coca3*cocc3 - unji*sibc*cocb3*cocc3)+pjmp*coca3*coca3

      ucor = xcod(zg,zi,0.0d0,uzj,uzzj)

      return
   end subroutine zcorp


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zfix here]
       subroutine zfixp(m,n,l,nx,np,n2,i,j,k,info,ilr,xa,xb,ya,yb,za,zb,t&
           ,hx,hy,hz,phi,x,y,z,zk,ui,cinfoz,cinfo,index2,fn,pjmp)


      implicit none
      _REAL_  x(0:m),y(0:n),z(0:l),phi(0:m,0:n,0:l) 
      integer index2(0:m,0:n,0:l) 
      _REAL_ ui(0:m,0:n,0:l),cinfoz(np*4*n2),fn(n2),cinfo(np*n2)
      integer m,n,l,nx,np,n2,i,j,k,info,ilr,nn,infodex
      _REAL_ xa,xb,ya,yb,za,zb,t,hx,hy,hz,xi,coca1,coca2,coca3
      _REAL_ cocb1,cocb2,cocb3,cocc1,cocc2,cocc3,curv,ftv,dft,ddft
      _REAL_ dsi,dsi2,hz1,hz2,sibb,sicc,sibc,zg,zk,yi,pjmp,ucor
      ! common /radius/r0

      hz2 = 2.0*hz
      hz1 = hz*hz
     infodex=3
      xi = x(i)
      yi = y(j)
      if(ilr == 1) then
         zg = z(k+1)
      else
         zg = z(k-1)
      end if

      coca1 = cinfoz(nx*np-np+4)
      coca2 = cinfoz(nx*np-np+5)
      coca3 = cinfoz(nx*np-np+6)
      cocb1 = cinfoz(nx*np-np+7)
      cocb2 = cinfoz(nx*np-np+8)
      cocb3 = cinfoz(nx*np-np+9)
      cocc1 = cinfoz(nx*np-np+10)
      cocc2 = cinfoz(nx*np-np+11)
      cocc3 = cinfoz(nx*np-np+12)
      sibb = cinfoz(nx*np-np+13)
      sicc = cinfoz(nx*np-np+14)
      sibc = cinfoz(nx*np-np+15)
      curv = cinfoz(nx*np-np+16)
      call pfjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,fn,x,y,z,ftv,dft,ddft,infodex,nx)
     nn=abs(index2(i,j,k))
     if(nn>0) then
     if(nint(cinfo(nn*np-np+17)) .eq. nint(cinfoz(nx*np-np+22))) then
     ftv=fn(nn)
     endif
     endif
      call zcorp(m,n,l,xa,xb,ya,yb,za,zb,t,hx,hy,hz,phi,x,y,z,zg &
            ,xi,yi,zk,coca1,coca2,coca3,cocb1,cocb2,cocb3 &
            ,cocc1,cocc2,cocc3,sibb,sicc,sibc,curv,ucor,ftv,dft,ddft,pjmp)

      dsi = -float(ilr)
      dsi2 = float(info)

      !       ui(i,j,k)=ui(i,j,k) + dsi2*rnu*0.5*(ucor+ucor2)/hz1
      ui(i,j,k)=ui(i,j,k) + dsi2*ucor

      return
   end subroutine zfixp

         subroutine interpp(m,n,l,np,n2,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,t,phi,fn,x,y,z &
         ,p,uvout,utout,usout,uout,pjmp,cinfo)
      
      implicit none
    
      integer nl,ne
      parameter( nl=27, ne=10 )

      !#-------------------------------------------------------------------------
      !# Subroutine interp3 interpolates u^+ of u_{ij} at (x1,y1) which is a point
      !# on the interface. The method is the weighted least squares interpolation.
      !#-- m,n:          The number of grid in the x- and y-direction.
      !#-- x1,y1:  A point on the interface where the interpolation take place.
      !#-- u(0:m,0:n):   The grid function to be interpolated.
      !#-------------------------------------------------------------------------

      

      _REAL_ b(nl),x(0:m),y(0:n),z(0:l),cinfo(np*n2),phi(0:m,0:n,0:l),p(0:m,0:n,0:l) 
       integer index2(0:m,0:n,0:l) 
       _REAL_ fn(n2),ss(nl),w2(nl,3),w3(nl),w3a(nl),w3b(nl) ,w3c(nl)
       integer ix(nl),iy(nl),iz(nl)
       integer m,n,l,np,n2,i,j,k,infodex,nx
       _REAL_ hx,hy,hz,xa,xb,ya,yb,za,zb,t
       _REAL_ uvout,usout,dft,ddft,utout,ftv,h,pjmp,sibb,sicc,sibc
       integer i1,i2,j1,k1,nn,nsub
       _REAL_ ts,unji,uout

      !# --- Get jump conditions --------------------------------------------
      infodex=4
      nx=abs(index2(i,j,k))
      call pfjump(m,n,l,n2,np,i,j,k,index2,hx,hy,hz &
         ,xa,xb,ya,yb,za,zb,fn,x,y,z,ftv,dft,ddft,infodex,nx)
     ! unji = ftv


      nn = abs(index2(i,j,k))
      unji = fn(nn)
     sibb=cinfo((nn-1)*np+13)
     sicc=cinfo((nn-1)*np+14)
     sibc=cinfo((nn-1)*np+15)
     nsub=numinterp(nn) 
      do i1=1,nsub
       ix(i1)=nint(interpinfo((nn-1)*npinterp+1,i1))
       iy(i1)=nint(interpinfo((nn-1)*npinterp+2,i1))
       iz(i1)=nint(interpinfo((nn-1)*npinterp+3,i1))
       ss(i1)=interpinfo((nn-1)*npinterp+4,i1)
       w2(i1,1)=interpinfo((nn-1)*npinterp+5,i1)
       w2(i1,2)=interpinfo((nn-1)*npinterp+6,i1)
       w2(i1,3)=interpinfo((nn-1)*npinterp+7,i1)
       w3(i1)=interpinfo((nn-1)*npinterp+8,i1)
       w3a(i1)=interpinfo((nn-1)*npinterp+9,i1)
       w3b(i1)=interpinfo((nn-1)*npinterp+10,i1)
       w3c(i1)=interpinfo((nn-1)*npinterp+11,i1)
      enddo
    
      !#-------- Begin to compute the correction term for u^{-} ------------

      uout=0.0d0
      uvout = 0.0d0
      utout = 0.0d0
      usout = 0.0d0




      do  i2=1, nsub
         i1 = ix(i2)
         j1 = iy(i2)
         k1 = iz(i2)
         ts = ss(i2)

         if (phi(i1,j1,k1)>0.0d0) then
           b(i2)=p(i1,j1,k1)*ts
            else
           if(solverorder .eq. 1) then 
           b(i2)=ts*(p(i1,j1,k1)+w2(i2,1)*unji)
           else
           b(i2)=ts*(p(i1,j1,k1)+w2(i2,1)*unji&
            +0.5*w2(i2,1)*w2(i2,1)*(unji*(sibb+sicc)+pjmp)& 
            -0.5*w2(i2,2)*w2(i2,2)*unji*sibb-0.5*w2(i2,3)*w2(i2,3)*unji*sicc&   
            +dft*w2(i2,1)*w2(i2,2)+ddft*w2(i2,1)*w2(i2,3)-unji*sibc*w2(i2,2)*w2(i2,3))
            endif
           endif  
         uout=uout+w3(i2)*b(i2)
         uvout= uvout + w3a(i2)*b(i2)
         utout= utout + w3b(i2)*b(i2)
         usout= usout + w3c(i2)*b(i2)

      enddo
      return
   end subroutine interpp
end module navier_stokes
