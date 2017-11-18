!######################################################################
! capmod: gateway for calling the routines that calculate the single-
!         particle basis representation of the CAP operator.
!######################################################################

module capmod

  implicit none  
  
  save
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter     :: dp=selected_real_kind(8)

contains
  
!######################################################################

  subroutine cap_mobas(gam,cap_mo)

    use channels
    use iomod
    use parameters
    use timingmod
    use monomial_analytic
    use import_gamess
    use misc, only: get_vdwr
    
    implicit none

    integer                               :: k
    real(dp)                              :: tw1,tw2,tc1,tc2
    real(dp), dimension(:,:), allocatable :: cap_mo
    type(gam_structure)                   :: gam
    
!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the MO representation of the &
         CAP operator'
    write(ilog,'(72a,/)') ('-',k=1,72)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(cap_mo(nbas,nbas))
    cap_mo=0.0d0

    allocate(vdwr(gam%natoms))
    vdwr=0.0d0

!----------------------------------------------------------------------
! Fetch the van der Waals radius array
!----------------------------------------------------------------------
    call get_vdwr(gam)
    
!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    if (icap.eq.1) then
       ! Monomial CAP, analytic evaluation of the CAP matrix elements
       call monomial_ana(gam,cap_mo,capord,capstr)
    else
       ! Numerical evaluation of the CAP matrix elements
       call numerical_cap(gam,cap_mo)
    endif

!----------------------------------------------------------------------
! Analysis of the MO representation of the CAP
!----------------------------------------------------------------------
    call analyse_cap_support(gam,cap_mo)

!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine cap_mobas

!######################################################################

  subroutine numerical_cap(gam,cap_mo)

    use channels
    use constants
    use parameters
    use basis_cap
    use misc, only: get_vdwr
    use import_gamess
    
    implicit none

    integer                                  :: nao,i,j,n,natom
    real(dp), dimension(nbas,nbas)           :: cap_mo
    real(dp), dimension(:,:), allocatable    :: cap_ao,smat,lmat
    real(dp), parameter                      :: ang2bohr=1.889725989d0
    real(dp), parameter                      :: dscale=3.5
    real(dp)                                 :: x,r
    complex(dp), dimension(:,:), allocatable :: cap_ao_cmplx
    type(gam_structure)                      :: gam
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nao=gam%nbasis
    natom=gam%natoms
        
    allocate(cap_ao(nao,nao))
    cap_ao=0.0d0

    allocate(cap_ao_cmplx(nao,nao))
    cap_ao_cmplx=czero

    allocate(smat(nao,nao))
    smat=0.0d0

    allocate(lmat(nao,nao))
    lmat=0.0d0
    
!----------------------------------------------------------------------
! Set the CAP type string used by the basis_cap module
!----------------------------------------------------------------------
    if (icap.eq.2) then
       ! Monomial CAP
       cap_type='monomial'
    else if (icap.eq.3) then
       ! Atom-centred monomial CAP
       cap_type='atom monomial'
    else if (icap.eq.4) then
       ! Moiseyev's non-local perfect CAP, single sphere
       cap_type='moiseyev'
    else if (icap.eq.5) then
       ! Moiseyev's non-local perfect CAP, atom-centered spheres
       cap_type='atom moiseyev'
    else if (icap.eq.6) then
       ! Cavity-like sigmoidal CAP
       cap_type='sigmoidal'
    endif

!----------------------------------------------------------------------
! CAP centre: for now we will take this as the geometric centre of
! the molecule
!----------------------------------------------------------------------
    cap_centre=0.0d0
    do n=1,natom
       do i=1,3
          cap_centre(i)=cap_centre(i) &
               +gam%atoms(n)%xyz(i)*ang2bohr/natom
       enddo
    enddo

!----------------------------------------------------------------------
! CAP starting radius: if a CAP box has not been specified by the user, 
! then we will take the start of the CAP to correspond to the greatest
! distance in any direction to the furthest most atom plus its van der
! Waals radius multiplied by dscale
!----------------------------------------------------------------------
    if (boxpar(1).eq.0.0d0) then
       ! The user has not specified a cap box
       cap_r0=-1.0d0
       do n=1,natom
          if (gam%atoms(n)%name.eq.'x') cycle
          do i=1,3
             x=gam%atoms(n)%xyz(i)*ang2bohr
             r=abs(x-cap_centre(i))+dscale*vdwr(n)
             if (r.gt.cap_r0) cap_r0=r
          enddo
       enddo
    else
       ! User specified box: take the maximum starting distance
       cap_r0=maxval(boxpar)
    endif

!----------------------------------------------------------------------
! CAP strength
!----------------------------------------------------------------------
    cap_strength=capstr

!----------------------------------------------------------------------
! CAP order
!----------------------------------------------------------------------
    cap_order=capord
    
!----------------------------------------------------------------------
! Calculate the AO representation of the CAP.
!
! Note that cap_evaluate calculates the AO representation of -iW, but
! we require the AO representation of W, hence the conversion after
! this subroutine is called.
!----------------------------------------------------------------------
    call cap_evaluate(gam,nrad(1),nang(1),nrad(2),nang(2),&
         cap_ao_cmplx,smat,lmat)

     do i=1,nao
        do j=1,nao
           cap_ao(i,j)=-aimag(cap_ao_cmplx(i,j))
        enddo
     enddo

!----------------------------------------------------------------------
! Transform the AO representation of the CAP to the MO representation
!----------------------------------------------------------------------
     cap_mo=matmul(transpose(ao2mo),matmul(cap_ao,ao2mo))
     
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(cap_ao)
    deallocate(cap_ao_cmplx)
    deallocate(smat)
    deallocate(lmat)
    
    return
    
  end subroutine numerical_cap

!######################################################################  

  subroutine analyse_cap_support(gam,cap_mo)

    use channels
    use constants
    use parameters
    use iomod
    use import_gamess
    use density
    
    implicit none

    integer                        :: c,i,j,p,q,unit
    integer, parameter             :: npnts=101
    real(dp), parameter            :: dr=0.25d0
    real(dp), dimension(nbas,nbas) :: cap_mo
    real(dp), dimension(3)         :: r
    real(dp), allocatable          :: aovalues(:)
    real(dp), allocatable          :: movalues(:)
    real(dp), allocatable          :: capval(:,:)
    type(gam_structure)            :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(aovalues(gam%nbasis))
    aovalues=0.0d0

    allocate(movalues(nbas))
    movalues=0.0d0

    allocate(capval(npnts,3))
    capval=0.0d0
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='mocapvals.dat',form='formatted',status='unknown')
    
!----------------------------------------------------------------------
! Calculate and output the values of the MO representation of the
! CAP at points along the x-, y-, and z-directions
!----------------------------------------------------------------------
    capval=0.0d0
    
    ! Loop over points
    do i=1,npnts
       
       ! Loop over directions
       do c=1,3
          
          ! Current coordinates
          r=0.0d0
          r(c)=(i-1)*dr
          
          ! Calculate the values of the AOs at the current point
          call get_ao_values(gam,aovalues,r)

          ! Calculate the values of the MOs at the current point
          movalues(1:nbas)=matmul(transpose(ao2mo),aovalues)
          
          ! Calculate the value of the MO representation of the
          ! CAP at the current point
          do p=1,nbas
             do q=1,nbas
                capval(i,c)=capval(i,c) &
                     +cap_mo(p,q)*movalues(p)*movalues(q)
             enddo
          enddo

       enddo

       ! Output the CAP values
       write(unit,'(F10.7,3(2x,ES15.8))') (i-1)*dr,(capval(i,j),j=1,3)
       
    enddo

!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(aovalues)
    deallocate(movalues)    
    
    return
    
  end subroutine analyse_cap_support
    
!######################################################################  
  
end module capmod
