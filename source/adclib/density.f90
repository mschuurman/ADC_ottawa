!######################################################################
! Routines for evaluating the electron density
!######################################################################

module density

  implicit none

  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter  :: dp=selected_real_kind(8)
  
contains

!######################################################################
! density_value: calcualates the value of the electron density at a
!                given position r for a given density matrix rho
!######################################################################
  
  function density_value(gam,rho,r) result(func)

    use parameters
    use import_gamess
    
    implicit none

    integer                        :: naos,p,q
    real(dp), dimension(nbas,nbas) :: rho
    real(dp), dimension(3)         :: r
    real(dp)                       :: func
    real(dp), allocatable          :: aovalues(:)
    real(dp), allocatable          :: movalues(:)
    type(gam_structure)            :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! No. AOs
    naos=gam%nbasis

    allocate(aovalues(naos))
    aovalues=0.0d0

    allocate(movalues(nbas))
    movalues=0.0d0
    
!----------------------------------------------------------------------
! Calculate the values of the AOs at the point r
!----------------------------------------------------------------------
    call get_ao_values(gam,aovalues,r,naos)

!----------------------------------------------------------------------
! Calculate the values of the MOs at the point r
!----------------------------------------------------------------------
    movalues(1:nbas)=matmul(transpose(ao2mo),aovalues)

!----------------------------------------------------------------------
! Calculate the value of the electron density
!----------------------------------------------------------------------
    func=0.0d0
    do p=1,nbas
       do q=1,nbas
          func=func+rho(p,q)*movalues(p)*movalues(q)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(aovalues)
    deallocate(movalues)
    
    return
    
  end function density_value

!######################################################################
! get_ao_values: calcuation of the values of all AOs at the point r
!######################################################################
  
  subroutine get_ao_values(gam,aovalues,r,dim)

    use parameters
    use import_gamess
    use gamess_internal

    implicit none

    integer                   :: dim,count,atm,sh,l,comp,pos,nx,ny,nz
    real(dp), dimension(dim)  :: aovalues
    real(dp), dimension(3)    :: r
    real(dp)                  :: x,y,z,angc
    real(dp), parameter       :: ang2bohr=1.889725989d0
    type(gam_structure)       :: gam
    
    ! AO counter
    count=0
    
    ! Loop over atoms
    do atm=1,gam%natoms

       ! Loop over shells
       do sh=1,gam%atoms(atm)%nshell

          ! Angular momentum quantum number, l
          l=gam%atoms(atm)%sh_l(sh)

          ! Loop over Cartesian components ({x,y,z},{xx,yy,...},etc)
          ! for the current l value
          do comp=1,gam_orbcnt(l)

             ! Increment the AO counter
             count=count+1

             ! Position of the current AO in
             ! the ang_n* arrays
             pos=ang_loc(l)+comp-1

             nx=ang_nx(pos)
             ny=ang_ny(pos)
             nz=ang_nz(pos)

             angc=ang_c(pos)

             ! Coordinate values relative to the atomic
             ! centres (in Bohr)
             x=r(1)-gam%atoms(atm)%xyz(1)*ang2bohr
             y=r(2)-gam%atoms(atm)%xyz(2)*ang2bohr
             z=r(3)-gam%atoms(atm)%xyz(3)*ang2bohr
                          
             ! AO value
             aovalues(count)=aoval(x,y,z,nx,ny,nz,angc,atm,sh,gam)
             
          enddo

       enddo

    enddo

    return
    
  end subroutine get_ao_values
    
!######################################################################

  function aoval(x,y,z,nx,ny,nz,angc,iatom,ishell,gam)

    use import_gamess
    use gamess_internal
    
    implicit none

    integer             :: nx,ny,nz,iatom,ishell,iprim,pindx1,pindx2
    real(dp)            :: aoval,x,y,z,zeta,coeff,rsq,angc
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Initialisation of the sum    
!----------------------------------------------------------------------
    aoval=0.0d0

!----------------------------------------------------------------------
! Distance squared
!----------------------------------------------------------------------
    rsq=x**2+y**2+z**2
    
!----------------------------------------------------------------------
! AO value
!----------------------------------------------------------------------
    ! Indices of the first and last primitives entering into the
    ! expansion of the current AO
    pindx1=gam%atoms(iatom)%sh_p(ishell)
    pindx2=gam%atoms(iatom)%sh_p(ishell+1)-1
    
    ! Loop over the primitives for the current AO
    do iprim=pindx1,pindx2

       ! Primitive exponent
       zeta=gam%atoms(iatom)%p_zet(iprim)

       ! Primitive coefficient
       coeff=gam%atoms(iatom)%p_c(iprim)
       
       ! Contribution of the current primitive to the AO
       aoval=aoval &
            + coeff &
            * angc &
            * (x**nx) &
            * (y**ny) &
            * (z**nz) &
            * exp(-zeta*rsq)
       
    enddo

    return
    
  end function aoval
  
!######################################################################
  
end module density
