! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
#include "pb_def.h"
#include "timer.h"

module poisson_boltzmann

   implicit none

#  include "pb_constants.h"

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = FOURPI*eps0*AMBER_ELECTROSTATIC2

   ! PBMD FD control variables

   logical :: outphi
   logical :: srsas
   logical :: scalerf
   logical :: outlvlset
   logical :: outmlvlset

   integer :: phiform 
   integer :: saopt
   integer :: sasopt
   integer :: dbfopt
   integer :: eneopt
   integer :: npbopt
   integer :: solvopt
   integer :: frcopt
   integer :: intopt
   integer :: bcopt
   integer :: smoothopt
   integer :: fold16
   integer :: isurfchg    
   integer :: membraneopt 
   integer :: poretype    
   integer :: augtoltype   
   integer :: rxm
   integer :: rym
   integer :: rzm
   integer :: xm
   integer :: ym
   integer :: zm
   integer :: xmym
   integer :: xmymzm
   integer :: nbuffer
   integer :: level
   integer :: nfocus
   integer :: fscale
   integer :: maxitn
   integer :: itn
   integer :: m, n
   integer :: savbcopt(MAXLEVEL)
   integer :: levelblock(MAXLEVEL)
   integer :: savxm(MAXLEVEL)
   integer :: savym(MAXLEVEL)
   integer :: savzm(MAXLEVEL)
   integer :: savxo(MAXLEVEL)
   integer :: savyo(MAXLEVEL)
   integer :: savzo(MAXLEVEL)
   integer :: savxmym(MAXLEVEL)
   integer :: savxmymzm(MAXLEVEL)
                          
   _REAL_ :: h
   _REAL_ :: rgox
   _REAL_ :: rgoy
   _REAL_ :: rgoz
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: fmiccg2  
   _REAL_ :: accept
   _REAL_ :: laccept
   _REAL_ :: wsor
   _REAL_ :: lwsor
   _REAL_ :: norm
   _REAL_ :: inorm
   _REAL_ :: xmax
   _REAL_ :: xmin
   _REAL_ :: ymax
   _REAL_ :: ymin
   _REAL_ :: zmax
   _REAL_ :: zmin
   _REAL_ :: gxmax
   _REAL_ :: gxmin
   _REAL_ :: gymax
   _REAL_ :: gymin
   _REAL_ :: gzmax
   _REAL_ :: gzmin
   _REAL_ :: savxbox(MAXLEVEL)
   _REAL_ :: savybox(MAXLEVEL)
   _REAL_ :: savzbox(MAXLEVEL)
   _REAL_ :: cxbox(MAXLEVEL)
   _REAL_ :: cybox(MAXLEVEL)
   _REAL_ :: czbox(MAXLEVEL)
   _REAL_ :: savh(MAXLEVEL)
   _REAL_ :: savgox(MAXLEVEL)
   _REAL_ :: savgoy(MAXLEVEL)
   _REAL_ :: savgoz(MAXLEVEL)
   _REAL_ :: offx
   _REAL_ :: offy
   _REAL_ :: offz
   _REAL_ :: fillratio

   _REAL_ :: epsin
   _REAL_ :: epsout
   _REAL_ :: epsmem
   _REAL_ :: epsmemb
   _REAL_ :: pbkappa
   _REAL_ :: istrng
   _REAL_ :: ivalence
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext

   _REAL_ :: mzmin
   _REAL_ :: mzmax
   _REAL_ :: mthick    
   _REAL_ :: mctrdz    
   _REAL_ :: poreradi  

   _REAL_ :: augctf    
   _REAL_ :: augtol    

   ! PBMD topology information

   integer              :: lastp
   integer              :: ngrdcrg

   integer, allocatable ::    icrd(:,:)
   integer, allocatable ::  grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)
 
   ! PBMD nblist information

   integer              :: maxnbr
   integer              :: maxnba
   _REAL_              :: cutres, cutnb, cutfd, cutsa
 
   integer, allocatable ::   nshrt(:)
   integer, allocatable ::     nex(:)
   integer, allocatable ::     iex(:,:)
   integer, allocatable :: iprshrt(:)
   integer, allocatable ::  iar1pb(:,:)
   _REAL_, allocatable :: cn1pb(:)
   _REAL_, allocatable :: cn2pb(:)
   _REAL_, allocatable :: cn3pb(:)

   ! PBMD cap water simulation information

   integer              :: mpopt
   integer              :: lmax
   integer              :: inatm
   integer              :: outwat
   integer              :: oution
   integer, allocatable :: outflag(:)
   integer, allocatable :: outflagorig(:)
   integer, allocatable :: mapout(:)
   integer, allocatable :: ibelly(:)
   _REAL_              :: sepbuf

   ! physical variables for energy and force calculations

   integer:: nbnd
   integer:: nbndx
   integer:: nbndy
   integer:: nbndz
   _REAL_, allocatable :: pos_crg(:,:,:)
   _REAL_, allocatable :: surf_crg(:,:)
   integer, allocatable :: ipos_crg(:,:,:)
   integer, allocatable :: crg_num(:)

   ! physical variable maps for numerical solutions

   _REAL_, allocatable ::     phi(:)
   _REAL_, allocatable ::      bv(:)
   _REAL_, allocatable ::   chgrd(:)
   _REAL_, allocatable ::    epsx(:)
   _REAL_, allocatable ::    epsy(:)
   _REAL_, allocatable ::    epsz(:)
   _REAL_, allocatable :: saltgrd(:)
   _REAL_, allocatable ::  ioncrg(:)

   ! geometry maps for dielectric interface
 
   integer, allocatable ::  insas(:)
   integer, allocatable :: atmsas(:)
   _REAL_, allocatable ::  lvlset(:)
   _REAL_, allocatable :: mlvlset(:)
   _REAL_, allocatable ::      zv(:)

   ! physical variable maps for FD force calculations

   _REAL_, allocatable ::     cphi(:)
   integer, allocatable ::  iepsav(:,:)
   integer, allocatable :: iepsavx(:,:)
   integer, allocatable :: iepsavz(:,:)
   integer, allocatable :: iepsavy(:,:)
   _REAL_, allocatable :: fedgex(:)
   _REAL_, allocatable :: fedgey(:)
   _REAL_, allocatable :: fedgez(:)

   ! saved phi array for pbmd
 
   _REAL_, allocatable :: xs(:)
   integer :: xsoffset

   ! ligand focusing options

   logical :: ligand
   character(len=256) ligandmask
   integer, allocatable :: liveflag(:)
   integer, allocatable :: realflag(:)
   integer :: ntrajmol
   _REAL_ :: buffer

   ! Multiple distributive fine grid geometry / Multiblock focusing

   logical :: multiblock                          !TRUE if specified multiblock
   logical :: firstleveldone                      !TRUE if level 1 is done once
   integer, allocatable :: blkxo(:)               !block origin in x dir
   integer, allocatable :: blkyo(:)               !block origin in y dir
   integer, allocatable :: blkzo(:)               !block origin in z dir
   integer, allocatable :: blkxlo(:)              !block lower bound in x
   integer, allocatable :: blkylo(:)              !block lower bound in y
   integer, allocatable :: blkzlo(:)              !block lower bound in z
   integer, allocatable :: blkxup(:)              !block upper bound in x
   integer, allocatable :: blkyup(:)              !block upper bound in y
   integer, allocatable :: blkzup(:)              !block upper bound in z
   integer, allocatable :: blknx(:)               !block grid length in x
   integer, allocatable :: blkny(:)               !block grid length in y
   integer, allocatable :: blknz(:)               !block grid length in z
   integer, allocatable :: blknxny(:)             !blknx . blkny
   integer, allocatable :: blknynz(:)             !blkny . blknz
   integer, allocatable :: blknxnz(:)             !blknx . blknz
   integer, allocatable :: blknxnynz(:)           !blknx . blkny . blknz
   integer              :: ngrdblkx               !# of grids per block in x
   integer              :: ngrdblky               !# of grids per block in y
   integer              :: ngrdblkz               !# of grids per block in z
   integer              :: xmblk                  !total blocks in x
   integer              :: ymblk                  !total blocks in y
   integer              :: zmblk                  !total blocks in z
   integer              :: ilower                 !
   integer              :: iupper                 !
   integer              :: jlower                 !
   integer              :: jupper                 !
   integer              :: klower                 !
   integer              :: kupper                 !
   integer              :: nfirstwpad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst0pad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst4s3f(MAXBLOCK+1) !Described in the text
   integer, allocatable :: lfocuswpad(:)          !Described in the text
   integer, allocatable :: lfocus0pad(:)          !Described in the text
   integer, allocatable :: lfocus4s3f(:)          !Described in the text
   integer, allocatable :: fineblockindex(:)      !for load balancing
   integer              :: iblkdebug              !
   _REAL_, allocatable :: blkgox(:)              !block origin in x
   _REAL_, allocatable :: blkgoy(:)              !block origin in y
   _REAL_, allocatable :: blkgoz(:)              !block origin in z
   _REAL_, allocatable :: coarsephi(:)           !saved 1st level solution

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of PBMD energy and forces
!     call pb_force( natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),&
!                    ix(i10),cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
subroutine pb_force( natom,nres,ntypes,npdec,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,eelrf )
    
   use solvent_accessibility, only : dprob, radi, radip, radip2, radip3, nzratm, &
#if defined SANDER || defined LIBPBSA
                                     sa_init, sa_driver, sa_free, sa_free_mb
#else
#ifdef MPI
                                     sa_init, sa_driver, sa_free, sa_free_mb,  &
                                     saslave_init
#else
                                     sa_init, sa_driver, sa_free, sa_free_mb
#endif /* MPI */
#endif /*SANDER or LIBPBSA*/
   use decomp, only : irespw, jgroup
   use pbtimer_module
   use random, only: amrset, amrand

   ! Common variables
    
#  include "../include/md.h"
#  include "pb_md.h"
#ifdef SANDER
#  include "../sander/box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#elif defined LIBPBSA
#  include "box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#else
#  include "box.h"
#ifdef MPI
   include "mpif.h"
#  include "parallel.h"
#else  /*MPI*/
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#endif /*MPI*/
#endif /* SANDER */
#  include "extra.h"
    
   ! Passed variables
    
   integer natom, nres, ntypes, npdec, ipres(*), iac(*), ico(*), natex(*)
   _REAL_ cn1(*), cn2(*), cg(natom), x(3,natom), f(3,natom)
   _REAL_ enb, eel, eelrf
 
   ! Local variables

   integer i, j, k
   integer iatm, proatm, atmfirst, atmlast
   integer atmind(natom)
   _REAL_ acg(natom)
   _REAL_ pbcutcap, pbxcap, pbycap, pbzcap
   _REAL_ eelrffd, eelrfmp
   _REAL_ pbfrc(3,natom)

   ! Local multiblock variables

   integer ipermute     !
   integer iblock       !
!  integer boolsorted
   integer mynblock     !
   integer lvl2begin    !
   integer orphanblkn   !
   integer guess_int    !
!#endif
   integer ierr         !
   integer taskpiece    !
   integer myblkstart   !
   integer myblkend     !
   integer tasknatom    !
   integer ihavedone    !
   ! readin info verification
   integer tmp_nsatm    !
   integer mingrdblk    !
   integer myldim       !
   logical indexmatched !

   ! This is not an efficient way to use the memory, need to estimate
   ! the size of the array.

   integer, allocatable :: tmpindex(:) !

   _REAL_ myh           !
   _REAL_ guess_float   !
   _REAL_ blk_eelrffd   ! temporary storage for eelrffd
   _REAL_ blk_eel       ! temporary storage for eel
   _REAL_ blk_enb       ! temporary storage for enb
#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
   _REAL_ , allocatable :: recvbuf1(:), recvbuf2(:), recvbuf3(:)
#endif /*MPI*/
#endif /*ndef SANDER or LIBPBSA*/

   ! end of multi-block

   logical localpbgrid

   ! Variables initialization
 
   enb = ZERO; eel = ZERO; eelrf = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO
   atmind = 0

   ! Variables initialization, multi-block

   myblkstart = 0
   myblkend   = 0
   ierr = 0
   acg = ZERO !making sure client(s) have clean acg

   firstleveldone = .false.

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
   call pbslave_init(natom)
   call saslave_init(natom)
#endif /* MPI */
#endif /*ndef SANDER or LIBPBSA*/

   ! End of multi-block initializatioin

   if ( ifcap /= 0 .and. (ifcap < 3 .or. ifcap > 5) ) then
      pbcutcap = cutcap+TWO; pbxcap = xcap; pbycap = ycap; pbzcap = zcap 
      radi(0) = pbcutcap; acrd(1,0) = pbxcap; acrd(2,0) = pbycap; acrd(3,0) = pbzcap
      radip3(1) = radi(0); nzratm(1) = 0
   else
      pbcutcap = ZERO; pbxcap = ZERO; pbycap = ZERO; pbzcap = ZERO
   end if

   ! part a.
   ! split atoms into internal/external and update nblist whenever a new pb grid
   ! is set up if requested
 
   call pbtimer_start(PBTIME_PBLIST)
   if ( pbgrid ) then
      if ( mpopt == 1 ) then
         ! multipole expansion
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,pbcutcap,sepbuf,x,ifcap)
      else if ( ifcap == 2 ) then 
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,cutcap,0.0d0,x,ifcap)
      else if ( ifcap == 5 ) then 
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart2(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                          inatm,outwat,oution,ipres,outflag, &
                          cutcap,x)
      else if ( ligand ) then
         ! This option will be visited if it is not the first time pb_force got
         ! called, however there would be something to do here once we've done
         ! MD.
         continue
      else
         ! Multiblock and other conditions go here
         outflag = 0
      end if
   else 
      if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) outflag = outflagorig
   end if

   if ( pqropt == 0 ) call pb_atmconv(mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg)

   ! part b. Set up nblist and grid
   ! This is for the global run for the coarse grid
   ! If ligand/multiple block is used, these will be updated later in docklist

   if ( ntnba == 1 .and. max(cutnb,cutsa,cutfd) > ZERO ) &
      call pb_atmlist(pbverbose,pbprint,pqropt,&
      maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,iar1pb,iprshrt,&
      cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,acrd(1,1))
   if ( ntnbr == 1 ) ntnbr = 0
   if ( ntnba == 1 ) ntnba = 0
   call pbtimer_stop(PBTIME_PBLIST)
 
   call pbtimer_start(PBTIME_PBSETUP)
   if ( mpopt /=2 .and. pbgrid ) then
      if ( ligand ) &
         call pb_atmpart3(pbverbose,pbprint,natom,buffer,xmin,xmax,ymin,ymax,&
              zmin,zmax,liveflag,realflag,outflag,x)
      if (ifcap == 2 .or. ifcap == 5) then
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,inatm,pbxcap,pbycap,pbzcap,pbcutcap)
      else
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,natom,pbxcap,pbycap,pbzcap,pbcutcap)
      end if
   end if
    
   call pbtimer_stop(PBTIME_PBSETUP)

   ! part c. compute grid-independent sas calculations for dielectric assignment
   ! when ifcap /= 0, no need to comptue sas
   call pbtimer_start(PBTIME_PBSAS)
   if ( srsas .and. ( ifcap == 0 .or. ifcap == 5 ) ) then
      if( ifcap == 5 ) then
         call sa_init(pbverbose,pbprint,natom,inatm,ifcap,dprob,radi,radip,radip2,outflag)
         ! the call here requires verification if we need to take care of multiblock as well
         call sa_driver(pbverbose,pbprint,pqropt,ipb,inp,natom,inatm,dosas,ndosas,npbstep,nsaslag,&
                        ligand, outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
      else
         call sa_init(pbverbose,pbprint,natom,natom,ifcap,dprob,radi,radip,radip2,outflag)
         call sa_driver(pbverbose,pbprint,pqropt,ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
                        (ligand .or. multiblock), outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
     end if
   end if

   call pbtimer_stop(PBTIME_PBSAS)

   call pbtimer_start(PBTIME_PBLIST)

   ! part d. Multi-block:
   ! Following part is for level 2 focusing task partition

   if ( multiblock ) then

      ! Current task partitioning: Predetermined distribution
      ! o o o o    For 4 simultaneous tasks doing 13 blocks, mytaskid
      ! o o o o    0, 1, 2 all get 3 blocks respectively. Mytaskid 3 
      ! o o o o    gets 4 blocks. Note: Queuing the blocks inside MPI
      !       o    is a pain in the neck, so it's not in our future
      !            plan.

      mynblock = levelblock(2)
      taskpiece  = levelblock(2)/numtasks
      orphanblkn = levelblock(2)-taskpiece*numtasks

      if ( numtasks < 2 ) then
      ! doing serial

         myblkstart = 1
         myblkend   = mynblock

      else if ( orphanblkn == 0 .or. mytaskid < numtasks-orphanblkn ) then
      ! doing parallel

         myblkstart =     mytaskid*taskpiece + 1
         myblkend   = (mytaskid+1)*taskpiece

      else if ( mytaskid == numtasks-orphanblkn ) then
      ! when orphanblkn == 1, it adds one block to the last task

         myblkstart =     mytaskid*taskpiece+1
         myblkend   = (mytaskid+1)*taskpiece+1

      else if ( mytaskid > numtasks-orphanblkn ) then

         myblkstart =     mytaskid*taskpiece+1+mytaskid+orphanblkn-numtasks
         myblkend   = (mytaskid+1)*taskpiece+1+mytaskid+orphanblkn-numtasks

      else

         write(6,'(a)') "exception caught, block division error for MPI"
         call mexit(6,1)

      end if

      if ( master ) then
         do i = 1, numtasks
            if ( orphanblkn == 0 .or. i < numtasks-orphanblkn ) then
               write(6,'(a,i6,a,i6,a)') "thread:",i,"is for",taskpiece  ,"fine blocks."
            else if ( i == numtasks-orphanblkn ) then
               write(6,'(a,i6,a,i6,a)') "thread:",i,"is for",taskpiece+1,"fine blocks."
            else if ( i > numtasks-orphanblkn ) then
               write(6,'(a,i6,a,i6,a)') "thread:",i,"is for",taskpiece+1,"fine blocks."
            end if
         end do
      end if
      taskpiece = myblkend - myblkstart + 1

      myh = savh(2)

      mingrdblk=min(ngrdblkx,ngrdblky,ngrdblkz)
      if ( buffer/myh > mingrdblk ) then
         write(6,'(a,i6,a,i6,a)') "Big padding",int(buffer/myh), &
                 "is unrealistic (the least grdblk is", &
                 mingrdblk,"), please give up."
      end if
      if (mingrdblk == 0) then
         write(*,*) "grdblkx or y, z contains zero."; stop
      end if

      ! Possible problem: mjhsieh forgot why there is a "+ 2"

      myldim=ceiling(2 + (buffer*2/myh/(mingrdblk-1)+3)**3) * natom
      myldim=min(myldim,mynblock*natom)

      ! Multithread situation is not considered so far, thus no benifit
      ! from multithread to reducing array sizes "for now".

      allocate(lfocuswpad(myldim),    stat=ierr); if (ierr /= 0) call mexit(6,1)
      allocate(lfocus0pad(2*natom+1), stat=ierr); if (ierr /= 0) call mexit(6,1)
      allocate(lfocus4s3f(myldim),    stat=ierr); if (ierr /= 0) call mexit(6,1)

      lfocus4s3f = 0 !initialization
      lfocuswpad = 0 !initialization
      lfocus0pad = 0 !initialization
      nfirst4s3f = 0 !initialization
      nfirstwpad = 0 !initialization
      nfirst0pad = 0 !initialization

      ! Following part initializes for level one

      nfirst4s3f(1) = 1
      nfirstwpad(1) = 1
      nfirst0pad(1) = 1
      ! Grid Partition for Level 2
      ! lfocuswpad, lfocus0pad, nfirstwpad and nfirst0pad got updated
      call pb_atmpart_mb(pbverbose, pbprint, natom, 1, mynblock, myh,&
              blknx,  blkny,  blknz,   blkxlo, blkylo, blkzlo,          &
              blkxup, blkyup, blkzup,  blkgox, blkgoy, blkgoz,          &
              x, lfocuswpad, lfocus0pad, lfocus4s3f,                    &
              nfirstwpad, nfirst0pad, nfirst4s3f, master)
      ! So the information needed for the solver are:
      !    myblkstart: the first block for focus finegrid in this thread
      !    myblkend:   the last block for focus finegrid in this thread
      !    lfocuswpad: serial number list of the atoms covered by the
      !                focus fine grid with the padding zone
      !    lfocus0pad: serial number list of the atoms covered by the
      !                focus fine grid without the padding zone
      !    lfocus4s3f: atom list for doing surface for blocks
      !    nfirstwpad: the first position in lfocuswpad that belongs to
      !                current block (level 2 blocks start from 2)
      !    nfirst0pad: the first position in lfocus0pad that belongs to
      !                current block (level 2 blocks start from 2)
      !    nfirst4s3f: the first position in lfocus4s3f

      tasknatom = nfirst0pad(mynblock+1)-nfirst0pad(1)
      if ( tasknatom .ne. natom ) then
         write(6, '(a,i6,a,i6)') 'pb_force(): Atom partition error',tasknatom,'/=',natom
         call mexit(6, 1)
      end if
   end if
   call pbtimer_stop(PBTIME_PBLIST)

   ! for focussing run, liveflag, outflag, realflag are updated
   ! atom list is updated next
   ! surface area is then updated
   if ( ligand ) then
      call pb_atmlist(pbverbose,pbprint,pqropt,maxnba,natom,ntypes,iac,ico,natex, &
              nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb, &
              cn2pb,cn3pb,cg,acrd)
      call sa_driver(pbverbose,pbprint,pqropt,ipb,inp,natom,natom,dosas,ndosas,    &
              npbstep,nsaslag,ligand,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
              nex,iex,.false.)

   ! BEGIN OF MULTIPLE BLOCK LOOP: LOOPING OVER BLOCKS

   else if ( multiblock ) then

      ! fineblockindex: storing the order of the multiblock for
      !                 load balancing or just sequential order

      allocate( fineblockindex(levelblock(2)), stat = ierr )
      do i = 1, levelblock(2)
         fineblockindex(i)=i
      end do

      ! This will shuffle the order, which is a fake load balancing.

      call amrset(8249562)
      if ( master ) then

         allocate(       tmpindex(levelblock(2)), stat = ierr )
         do i = 1, levelblock(2)
            call amrand(guess_float) ! from 0 to 1
            guess_float = guess_float * (levelblock(2)+1-i)
            guess_int   = ceiling(guess_float)
            ipermute    = fineblockindex(guess_int)
            fineblockindex(guess_int) = fineblockindex(levelblock(2)+1-i)
            fineblockindex(levelblock(2)+1-i) = ipermute
         end do
         tmpindex=fineblockindex

         ! distributing blocks among available nodes

         k = 1
         do i = 0, numtasks-1
            do j = 1, levelblock(2)
               if ( mod(j,numtasks) == i ) then
                  fineblockindex(k)=tmpindex(j)
                  k = k + 1
               end if
            end do
         end do
         deallocate( tmpindex, stat = ierr )
         ihavedone = 0
      end if

      ! broadcast the ownership of blocks to each thread

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
      call MPI_BCAST(fineblockindex(1),levelblock(2),MPI_INTEGER,0,CommSANDER,ierr)
      REQUIRE(ierr==0)
#endif /*MPI*/
#endif /*ndef SANDER or LIBPBSA*/

      ! add FD reaction field energy/force

      call pbtimer_start(PBTIME_PBFDFRC)
      blk_eelrffd = ZERO
      blk_eel     = ZERO
      blk_enb     = ZERO
      allocate(coarsephi(savxmymzm(1)), stat = ierr)
      blkloop: do iblock=1, levelblock(2)
         iblkdebug = iblock
         indexmatched = .false.
         matchingindex: do i = myblkstart, myblkend
            if ( iblock == fineblockindex(i) ) then
               indexmatched = .true.
               cycle matchingindex
            end if
         end do matchingindex
         if ( .not. indexmatched ) cycle blkloop
         liveflag = 0 ! 1 if inside the block  w/o pad
         realflag = 0 ! 1 if inside the block with pad
         outflag  = 1 ! 1 if outside the geometry for surface
         do i = nfirst0pad(iblock), nfirst0pad(iblock+1)-1
            liveflag( lfocus0pad(i) ) = 1
         end do
         do i = nfirstwpad(iblock), nfirstwpad(iblock+1)-1
            realflag( lfocuswpad(i) ) = 1
         end do
         do i = nfirst4s3f(iblock), nfirst4s3f(iblock+1)-1
            outflag ( lfocus4s3f(i) ) = 0
         end do
         savxm(2)=blknx(iblock)
         savym(2)=blkny(iblock)
         savzm(2)=blknz(iblock)
         savxmym(2)=blknxny(iblock)
         savxmymzm(2)=blknxnynz(iblock)
         savgox(2)=blkgox(iblock)
         savgoy(2)=blkgoy(iblock)
         savgoz(2)=blkgoz(iblock)
         call pb_atmlist(pbverbose,pbprint,pqropt,maxnba,natom,ntypes,iac,ico,natex, &
                 nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb, &
                 cn2pb,cn3pb,cg,acrd)

         ! In this implementation of multiblock focusing, each block is treated
         ! as exactly like a ligand box, so we hijack the ligand flag for now.
         ! Here is the interface before the hijack:
         ! call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,  &
         !       npbstep,nsaslag,ligand,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
         !       nex,iex)
         call sa_driver(pbverbose,pbprint,pqropt,ipb,inp,natom,natom,dosas,ndosas,    &
                 npbstep,nsaslag,.true.,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
                 nex,iex,.false.)

         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,    &
                 npdec,idecomp,irespw,ipres,jgroup,ibgwat,ibgion,pbfrc,        &
                 enb,eelrffd,npbstep,npbgrid,nstlim)

         blk_eelrffd = blk_eelrffd + eelrffd
         ihavedone = ihavedone + 1
         if ( srsas .and. ihavedone < taskpiece ) then
            call sa_free_mb( dosas,ndosas )
         end if

         ! cutnb > 0 is required in multiblock

         if( ifcap == 2 .or. ifcap == 5) then
            call pb_directwtcut(natom,inatm,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                    pbfrc,eel,enb)
         else
            call pb_directwtcut(natom,natom,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                    pbfrc,eel,enb)
         end if 
         blk_eel = blk_eel + eel
         blk_enb = blk_enb + enb
      end do blkloop
      deallocate (coarsephi, stat = ierr)
      outflag = 0
      realflag = 1
      eelrffd = blk_eelrffd
      eel = blk_eel
      enb = blk_enb
      ! PHI MAP SHOULD BE EXPORTED IN FDFRC, NOT HERE.
      ! To do so xs is required for assembling the phi
      call pbtimer_stop(PBTIME_PBFDFRC)

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
      allocate(recvbuf1(numtasks), stat = ierr); recvbuf1 = ZERO
      allocate(recvbuf2(numtasks), stat = ierr); recvbuf2 = ZERO
      allocate(recvbuf3(numtasks), stat = ierr); recvbuf3 = ZERO
      call MPI_GATHER(eelrffd,  1, MPI_DOUBLE_PRECISION, recvbuf1, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      call MPI_GATHER(eel    ,  1, MPI_DOUBLE_PRECISION, recvbuf2, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      call MPI_GATHER(enb    ,  1, MPI_DOUBLE_PRECISION, recvbuf3, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      if (master) then
         eelrffd = SUM(recvbuf1)
         eel     = SUM(recvbuf2)
         enb     = SUM(recvbuf3)
      end if
#endif /* MPI */
#endif /*ndef SANDER or LIBPBSA*/

   end if

   ! HERE GOES THE ENERGY SUMMATION (parallel)
   ! END OF MULTIPLE BLOCK LOOP

   ! part e. add FD reaction field energy/force

   call pbtimer_start(PBTIME_PBFDFRC)
   if ( multiblock ) then
      ! For multiblock focusing, FD force is already done in the previous
      ! section. But this is not the case for ligand.
      continue
   else if ( epsout /= epsin .and. mpopt /= 2 ) then
      ! In the case of ifcap == 2,5, only map crg within cap to grid (atmlast == inatm),
      ! else map all crg (atmlast == natom)
      if( ifcap == 2 .or. ifcap == 5 ) then
         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,inatm,npdec,idecomp,irespw, &
                 ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
      else
         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                 ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim)
      end if
   else if (epsout == epsin) then
      ! Apparently this is for expert use only
         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                 ipres,jgroup,ibgwat,ibgion,pbfrc,enb,eelrffd,npbstep,npbgrid,nstlim) 
   else
      write(6,'(a)') 'PB Bomb in pb_force(): Unknown FD force call option'
      call mexit(6,1)
   end if
   call pbtimer_stop(PBTIME_PBFDFRC)

   ! clean up for sas calculations

   call pbtimer_start(PBTIME_PBSETUP)
   if ( srsas .and. (ifcap == 0 .or. ifcap == 5) ) then
      call sa_free( dosas,ndosas,.false. )
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   if ( saopt < 0 ) return

   ! part f. add MP reaction field energy/forces when ifcap /= 0
 
   call pbtimer_start(PBTIME_PBMP)
   if ( mpopt /= 0 .and. epsout /= epsin ) then
      if ( mpopt == 1 ) then      ! multipole expansion for boundary atoms
         atmfirst = inatm + 1
         atmlast  = natom
      else if ( mpopt == 2 ) then ! multipole expansion for all atoms
         atmfirst = 1
         atmlast  = natom
      end if
      call pb_mpfrc(natom,atmfirst,atmlast,lmax,pbcutcap,pbxcap,pbycap,pbzcap,&
              epsin,epsout,acrg,acrd(1,1),pbfrc,eelrfmp)
   end if
   call pbtimer_stop(PBTIME_PBMP)
 
   ! part g. add direct coulombic and nonbonded forces
    
   call pbtimer_start(PBTIME_PBDIRECT)
   if ( multiblock .or. eneopt == 4) then
      ! For multiblock focusing or new p3m, pb_direct SHOULD BE DONE inside the loop
      continue
   else if ( pqropt == 0 .and. cutnb == ZERO ) then
      call pb_directnocut(natom,proatm,inatm,ipres,ibgwat,ienwat,ibgion,ienion,ntypes,eneopt,idecomp,ifcap, &
              iac,ico,nex,iex,cn1,cn2,acg,acrd(1,1),pbfrc,eel,enb)
   else
      if( ifcap == 2 .or. ifcap == 5) then
         call pb_directwtcut(natom,inatm,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                 pbfrc,eel,enb)
      else if ( pqropt == 0 ) then
         call pb_directwtcut(natom,natom,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                 pbfrc,eel,enb)
      end if
   end if
   call pbtimer_stop(PBTIME_PBDIRECT)

   ! part h. returning:
   ! i. returning energies

   eel = eel * eps0/epsin
   if ( eneopt == 1 .and. (bcopt < 6 .or. bcopt > 9) ) then 
      eel = eel + eelrffd + eelrfmp
      eelrf = ZERO
   else if ( eneopt == 4 ) then
      eel = eelrffd
      eelrf = ZERO
   else
      eelrf = eelrffd + eelrfmp
   end if

   ! ii. returning forces
   ! only do this for normal amber file input

   if ( pqropt == 0 ) then
      if ( ifcap == 2 .or. ifcap == 5 ) then
         atmlast = inatm
      else
         atmlast = natom
      end if

      do iatm = 1, natom
         f(1,iatm) = f(1,iatm) + pbfrc(1,mapout(iatm))
         f(2,iatm) = f(2,iatm) + pbfrc(2,mapout(iatm))
         f(3,iatm) = f(3,iatm) + pbfrc(3,mapout(iatm))
      end do
   end if

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ convert passed coordinates and charges to the internal format
subroutine pb_atmconv( mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
   
   integer mpopt, ifcap, natom, ibgwat, ienwat, ibgion, ienion
   integer atmind(natom), ipres(*) 
   _REAL_ x(3,natom), cg(natom), acg(natom)
    
   ! Local variables
    
   integer i, j, ifirst, ilast, iatm, ires, num
    
   if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) then
       
      ! copy reordered coord/charge to private arrays for pb/mp or ifcap == 2,5
      ! protein atoms go into internal portion (so do IONS!)
       
      ifirst = 1; ilast = ipres(ibgwat)-1

      do iatm = ifirst, ilast
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm); atmind(iatm) = iatm
         mapout(iatm) = iatm
      end do
       
      ! water atoms go into internal/external portion
       
      ifirst = ipres(ibgwat); ilast = natom

      i = ifirst; j =  inatm + 1
      do iatm = ifirst, ilast
         if ( outflag(iatm) == 0 ) then
            acrd(1,i   ) = x(1,iatm); acrd(2,i   ) = x(2,iatm); acrd(3,i   ) = x(3,iatm)
            acrg(i   ) = cg(iatm)/18.2223d0; acg(i   ) = cg(iatm); atmind(i   ) = iatm
            mapout(iatm) = i
            i = i + 1
         else
            acrd(1,j   ) = x(1,iatm); acrd(2,j   ) = x(2,iatm); acrd(3,j   ) = x(3,iatm);
            acrg(j   ) = cg(iatm)/18.2223d0; acg(j   ) = cg(iatm); atmind(j   ) = iatm
            mapout(iatm) = j
            j = j + 1
         end if
      end do

      ! store original outflag array and prepare an updated one for water atoms

      outflagorig = outflag
      do iatm = ifirst, ilast
         if( iatm <= inatm ) then
            outflag(iatm) = 0
         else
            outflag(iatm) = 1
         end if
      end do
       
   else
       
      do iatm = 1, natom
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
         mapout(iatm) = iatm
      end do
   end if
 
end subroutine pb_atmconv

end subroutine pb_force
#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ MPI slave initiation
!  Mengjuei Hsieh, University of California Irvine
subroutine pbslave_init(natom)
   implicit none
#  include "flocntrl.h"
#  include "pb_md.h"
   include "mpif.h"
#  include "parallel.h"
#  include "extra.h"
   integer natom
   _REAL_ green(0:40, 0:40, 0:40)
   common /blk_green/ green
   !LOCAL
   integer ierr
   !MPI initialization VERY BAD IMPLEMENTATION
   if ( .not. master ) then
      pbprint = .false.
      allocate(liveflag(    natom),STAT=ierr)
      allocate(realflag(    natom),STAT=ierr)
      allocate(outflag (    natom),STAT=ierr)
      allocate(outflagorig( natom),STAT=ierr)
      allocate(  iar1pb(4,0:natom),STAT=ierr)
      allocate(   nshrt(  0:natom),STAT=ierr)
      allocate(    acrd(3,  natom),STAT=ierr)
      allocate(    gcrd(3,  natom),STAT=ierr)
      allocate(    icrd(3,  natom),STAT=ierr)
      allocate(    acrg(    natom),STAT=ierr)
      allocate(    gcrg(8,  natom),STAT=ierr)
      allocate(  grdcrg(3,8*natom),STAT=ierr)
      allocate( qgrdcrg(  8*natom),STAT=ierr)
      allocate(  mapout(    natom),STAT=ierr)
      allocate(     nex(    natom),STAT=ierr)
      allocate(     iex(64, natom),STAT=ierr)
   endif
   call MPI_BCAST(    ligand,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(multiblock,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     srsas,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxnbr,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxnba,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    nfocus,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    fscale,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   nbuffer,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   solvopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     bcopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    eneopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    intopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxitn,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblkx,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblky,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblkz,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  nshrt(0),    natom+1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    nex(1),      natom,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  iex(1,1),      natom,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    do_dir,BC_FLOCNTRL,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    buffer,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    cutres,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutnb,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutsa,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutfd,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offx,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offy,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offz,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST( fillratio,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    fmiccg,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    accept,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    pbtemp,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ivalence,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      wsor,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     lwsor,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     epsin,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    epsout,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST( acrd(1,1),natom*3,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(green(0,0,0), 9261,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
   call MPI_BCAST(savbcopt(1),nfocus,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   savh(1), nfocus,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
   if ( .not. master ) then
      allocate(iprshrt (maxnba),stat=ierr)
      allocate(cn1pb   (maxnba),stat=ierr)
      allocate(cn2pb   (maxnba),stat=ierr)
      allocate(cn3pb   (maxnba),stat=ierr)
   end if
end subroutine pbslave_init
#endif /*def MPI*/
#endif /*ndef SANDER or LIBPBSA*/

end module poisson_boltzmann
