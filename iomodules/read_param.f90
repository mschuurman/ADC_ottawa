module read_param

  use constants
  use parameters
  use misc
  use channels

  implicit none
  
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine load_gamess(chkpt_file,log_file)
  use accuracy
  use parameters
  use import_gamess
  use math
  character(len=72),intent(inout)   :: chkpt_file,log_file
  type(gam_structure)               :: gamess_info       ! Default GAMESS
  integer                           :: j,naos
  real(xrk)                         :: q

  ! determine the number of aos, mos, and read in the mos
  write(ilog,'(a)') 'chk='//trim(chkpt_file)
  call accuracyInitialize
  call gamess_load_orbitals(file=trim(chkpt_file),structure=gamess_info)
  write (ilog,"(/'Loaded GAMESS checkpoint file ',a/)") trim(chkpt_file)

  nBas  = gamess_info%nvectors
  naos  = gamess_info%nbasis

  allocate(orbSym(nBas),e(nBas))

  ! determine various electronic structure variables
  call read_gamess_output(nBas,nelec,nCen,nIrr,e,orbSym,labSym,pntgroup)
  write (ilog,"(/'Loaded GAMESS log file ',a/)") trim(log_file)

  ! load MO integrals into memory
  call load_mo_integrals(gamess_info)
  write (ilog,"(/'Loaded MO integrals into memory ',a/)")

  write(ilog,102) "HF energy", Ehf, "a.u."
  write(ilog,100) 'nr.','sym','sym label','occ','orbital energy'
  write(ilog,'(10x,50("-"))')
  do j=1,nBas
     write(ilog,101) j,orbSym(j),labSym(orbSym(j)),occNum(j),e(j)
  end do  

100 FORMAT(/,10x,A3,5x,A3,5x,A9,5x,A3,5x,A16)
101 FORMAT(/,10x,I3,5x,I3,5x,A3,5x,F3.0,5x,F16.10)
102 FORMAT(/,3("-"),A30,5x,F16.10,1x,A4)
end subroutine load_gamess
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine read_user()

  integer :: k

  NAMELIST /USER/ debug,nirrep,hcentre,dlim,minc,stiprilev,method,matvec,idiag,fdiag,fmethod,WHAT, &
       davname,lancname,mspacewi,mspacewf,davstates,numinista,chrun,chrun2,NSYMA,ELECTRIC_FIELD,POLARIZATION, &
       eupper,elower,readband,tranmom,norder,info,ninista,statenumber,nirrep2,tranmom2,denord,GO,DIPOLESYM,&
       lcvs,icore,lfakeip,ifakeorb,expfakeip,moType,ld2,lcvsfinal,lscf,dmatmem

  NAMELIST /LNZLST/ ncycles,maxmem,memx,mode,nprint,maxiter,wthr,erange,unit,fparm,lmain,dmain,davtol,&
       ladc1guess

  hcentre(:)=-1

  dlim=250.0d0

  moType = 'incore'

  ninista=-1
  debug=.false.

  ladc1guess=.false.
  davtol=1d-7

  lcvs=.false.  
  lcvsfinal=.false.
  icore=0

  lfakeip=.false.
  ifakeorb=0

  ld2=.false.

  lscf=.true.

  dmatmem=250.0d0

  READ(*,USER)
  READ(*,LNZLST)

  stvc_lbl(:)=-1

  iparm(1)=maxiter
  hinit=hcentre(1)
  lancstates=ncycles*lmain
  
  ncore=0
  do k=1,nhcentre
     if (icore(k).gt.0) ncore=ncore+1
  enddo

  if (lfakeip.and.ifakeorb.eq.0) then
     write(ilog,'(/,2x,a,/)') 'Index of the fake IP diffuse orbital has &
          not been given'
     STOP
  endif

  write(ilog,*) 'Orbitals on the central atom', hcentre(0),':',hcentre(1:hcentre(0))
  
  write(ilog,*) "User data read"
  
end subroutine read_user


!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine check_user()
  

  if ((nirrep .lt. 1) .or. (nirrep .gt. 8)) then
     write(ilog,*) "IRREP IS NOT SET"
     stop
!  elseif ((hinit .le. 0) .or. (hinit .gt. nBas)) then
!     write(ilog,*) "INITIAL HOLE ORBITAL NUMBER IS NOT SET"
!     stop
!  elseif ((hinit .le. 0) .or. (hinit .gt. nBas)) then   
!     write(ilog,*) "INITIAL PARTICLE ORBITAL NUMBER IS NOT SET"
!     stop 
!  elseif ((method .lt. 0) .or. (method .gt. 20)) then
!  elseif (method .gt. 20) then
!     write(ilog,*) "COMPUTATION METHOD IS NOT CHOSEN"
!     stop
!  elseif ((numinista .lt. 0) .or. (numinista .gt. davstates)) then
!     write(ilog,*) "THE NUMBER OF THE INITIAL STATES IS WRONG"
!     stop   
  end if
     


end subroutine check_user

!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine rearrange_occ()

!!$ The subroutine rearranges MO's such that occupied MO's precede unoccupied.
!!$ roccNum contains old indices, nOcc number of occup. orbs.
  
  integer :: i,j,k
  integer, dimension(nBas):: indx
  
100 FORMAT(/,3("-"),A50,3x,I4)
101 FORMAT(/,60("-"),/)

103 FORMAT(/,10x,A3,5x,A3,5x,A3,5x,A16)
104 FORMAT(/,10x,I3,5x,I3,5x,F3.1,5x,F16.10)
  
  j=0
  k=0
  do i=1,nBas
     if (occNum(i).eq.2) then
        j=j+1
!        roccNum(j)=i
     else
        k=k+1
!        roccNum(nAct-k+1)=i
     end if
  end do
  
  nOcc=j
  nVirt=k
  
  write(ilog,'(/,60("*"),/)')
  if ((k+j) .ne. nBas) then
     write(ilog,*) "(rearrange_occ):Orbital number is wrong. Bailing out ..."
     stop
  end if
  
  write(ilog,100) "hole space size",nOcc
  write(ilog,100) "particle space size",nVirt 
  write(ilog,101)
  
!!$  call dsortqx("A",nBas,e(:),1,indx(:))
  
!  call dsortindxa1('A',nBas,e(:),indx(:))
  indx = (/(i, i=1,nBas )/)
  write(ilog,*)'index=',indx
  roccNum(:)=indx(:)
  
end subroutine rearrange_occ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rdorbsym

  implicit none
  
  integer            :: ilog
  character(len=120) :: string

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
  ilog=735
  open(ilog,file='gamess.log',form='formatted',status='old')

!-----------------------------------------------------------------------
! Read to the MO section
!-----------------------------------------------------------------------
10 read(ilog,'(a)') string
  if (index(string,'EIGENVECTORS').eq.0) goto 10
  read(ilog,*)
  read(ilog,*)

!-----------------------------------------------------------------------
! Read the MO symmetries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
  close(ilog)


  return

end subroutine rdorbsym
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module read_param
