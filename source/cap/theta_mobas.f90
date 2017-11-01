!######################################################################
! thetamod: calculation of the MO representation of the projector,
!           Theta, onto the CAP region (also known as the
!           characteristic function)
!######################################################################

module thetamod

  implicit none

  save
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter     :: dp=selected_real_kind(8)

contains
  
!######################################################################

  subroutine theta_mobas(gam,theta_mo)

    use channels
    use iomod
    use parameters
    use timingmod
    use monomial_analytic
    use import_gamess

    implicit none

    integer                               :: k
    real(dp)                              :: tw1,tw2,tc1,tc2
    real(dp), dimension(:,:), allocatable :: theta_mo
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the MO representation of the &
         projector onto the CAP region'
    write(ilog,'(72a,/)') ('-',k=1,72)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(theta_mo(nbas,nbas))
    theta_mo=0.0d0

!----------------------------------------------------------------------    
! Calculate the MO representation of the projector onto the CAP region
!----------------------------------------------------------------------
    if (icap.eq.1) then
       call monomial_ana(gam,theta_mo,0,1.0d0)
    else
       call numerical_theta(gam,theta_mo)
    endif
    
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"

    return
    
  end subroutine theta_mobas

!######################################################################

  subroutine numerical_theta(gam,theta_mo)

    use channels
    use constants
    use parameters
    use basis_cap
    use misc, only: get_vdwr
    use iomod
    use import_gamess
    
    implicit none

    integer                                  :: nao,i,j,n,natom
    real(dp), dimension(nbas,nbas)           :: theta_mo
    real(dp), dimension(:,:), allocatable    :: theta_ao,smat,lmat
    real(dp), parameter                      :: ang2bohr=1.889725989d0
    real(dp), parameter                      :: dscale=3.5
    real(dp)                                 :: x,r
    complex(dp), dimension(:,:), allocatable :: theta_ao_cmplx
    type(gam_structure)                      :: gam
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nao=gam%nbasis
    natom=gam%natoms
    
    if (.not.allocated(vdwr)) then
       allocate(vdwr(natom))
       vdwr=0.0d0
    endif
       
    allocate(theta_ao(nao,nao))
    theta_ao=0.0d0

    allocate(theta_ao_cmplx(nao,nao))
    theta_ao_cmplx=czero

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
    else if (icap.eq.6) then
       ! Cavity-type sigmoidal CAP
       cap_type='sigmoidal'       
    else
       errmsg='Flux analysis is not yet supported for ECS-type &
            absorbing potentials'
       call error_control
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
       ! The user has not specified a cap box, 
       call get_vdwr(gam)
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
! CAP strength: 1 for the calculation of the MO representation of
! Theta
!----------------------------------------------------------------------
    cap_strength=1.0d0

!----------------------------------------------------------------------
! CAP order: 0 for the calculation of the MO representation of
! Theta
!----------------------------------------------------------------------
    cap_order=0
    
!----------------------------------------------------------------------
! Calculate the AO representation of Theta.
!
! Note that cap_evaluate calculates the AO representation of -iTheta,
! but we require the AO representation of Theta, hence the conversion
! after this subroutine is called.
!----------------------------------------------------------------------
    call cap_evaluate(gam,nrad(1),nang(1),nrad(2),nang(2),&
         theta_ao_cmplx,smat,lmat)

     do i=1,nao
        do j=1,nao
           theta_ao(i,j)=-aimag(theta_ao_cmplx(i,j))
        enddo
     enddo

!----------------------------------------------------------------------
! Transform the AO representation of Theta to the MO representation
!----------------------------------------------------------------------
     theta_mo=matmul(transpose(ao2mo),matmul(theta_ao,ao2mo))
     
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(theta_ao)
    deallocate(theta_ao_cmplx)
    deallocate(smat)
    deallocate(lmat)
    
    return
    
  end subroutine numerical_theta

!######################################################################

end module thetamod
