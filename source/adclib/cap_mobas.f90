!######################################################################
! capmod: routines for the calculation of the MO representation of
!         the CAP operator. Employs numerical integration performed
!         using Becke's partitioning scheme (JCP 88, 2547 (1988)).
!         Makes use of the NUMGRID libraries of Radovan Bast.
!######################################################################

module capmod

  use constants
  use numgrid
  use, intrinsic :: iso_c_binding, only: c_ptr
  
  implicit none  
  
  save
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter    :: dp=selected_real_kind(8)

  ! Conversion factors
  real(dp), parameter   :: ang2bohr=1.889725989d0
  
  ! NUMGRID arrays and variables
  integer               :: min_num_angular_points
  integer               :: max_num_angular_points
  integer               :: num_points
  integer               :: num_centers
  integer               :: num_outer_centers
  integer               :: num_shells
  integer               :: num_primitives
  integer, allocatable  :: center_elements(:)
  integer, allocatable  :: shell_centers(:)
  integer, allocatable  :: shell_l_quantum_numbers(:)
  integer, allocatable  :: shell_num_primitives(:)
  integer, allocatable  :: outer_center_elements(:)
  real(dp), allocatable :: outer_center_coordinates(:)
  real(dp), allocatable :: primitive_exponents(:)
  real(dp)              :: radial_precision
  real(dp), allocatable :: center_coordinates(:)
  real(dp), pointer     :: grid(:)
  type(c_ptr)           :: context
  
contains
  
!######################################################################

  subroutine cap_mobas(gam)

    use channels
    use constants
    use parameters
    use timingmod
    use import_gamess
    
    implicit none

    integer             :: k
    real(dp)            :: tw1,tw2,tc1,tc2
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('-',k=1,72)
    write(ilog,'(2x,a)') 'Calculating the MO representation of the &
         CAP operator'
    write(ilog,'(72a)') ('-',k=1,72)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Set up the integration grid
!----------------------------------------------------------------------
    call initialise_intgrid(gam)

!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    call mo_cap_matrix(gam)

!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
!----------------------------------------------------------------------
! Deallocate arrays and destroy the grid context
!----------------------------------------------------------------------
    call finalise_intgrid
    
    return
    
  end subroutine cap_mobas

!######################################################################

  subroutine initialise_intgrid(gam)

    use channels
    use constants
    use parameters
    use import_gamess
    
    implicit none

    integer             :: i,j,k,l
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Grid parameters
!----------------------------------------------------------------------
    ! Grid parameters
    radial_precision=1.0d-14
    min_num_angular_points=170
    max_num_angular_points=350

!----------------------------------------------------------------------
! Atomic coordinates (in Bohr). Note that the gam_structure derived
! type holds the coordinates in Angstrom.
!----------------------------------------------------------------------
    ! Number of atoms
    num_centers=gam%natoms

    ! Atom coordinates
    allocate(center_coordinates(num_centers*3))

    k=0
    do i=1,gam%natoms
       do j=1,3
          k=k+1
          center_coordinates(k)=gam%atoms(i)%xyz(j)*ang2bohr
       enddo
    enddo

!----------------------------------------------------------------------
! Atom types
!----------------------------------------------------------------------
    allocate(center_elements(num_centers))
    do i=1,gam%natoms
       center_elements(i)=int(gam%atoms(i)%znuc)
    enddo

!----------------------------------------------------------------------
! Basis information
!----------------------------------------------------------------------
    num_outer_centers=0

    ! Number of shells
    num_shells=0
    do i=1,gam%natoms
       num_shells=num_shells+gam%atoms(i)%nshell
    enddo

    ! Shell centers
    allocate(shell_centers(num_shells))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_centers(k)=i
       enddo
    enddo

    ! Angular momentum quantum numbers for each shell
    allocate(shell_l_quantum_numbers(num_shells))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_l_quantum_numbers(k)=gam%atoms(i)%sh_l(j)
       enddo
    enddo

    ! Number of primitives for each shell
    allocate(shell_num_primitives(num_shells))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_num_primitives(k)=gam%atoms(i)%sh_p(j+1)-gam%atoms(i)%sh_p(j)
       enddo
    enddo

    ! Total number of primitives
    num_primitives=sum(shell_num_primitives)

    ! Primitive exponents
    allocate(primitive_exponents(num_primitives))
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
             k=k+1
             primitive_exponents(k)=gam%atoms(i)%p_zet(l)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Create the grid
!----------------------------------------------------------------------
    context=numgrid_new_context()

    call numgrid_generate_grid(context,&
         radial_precision,&
         min_num_angular_points,&
         max_num_angular_points,&
         num_centers,&
         center_coordinates,&
         center_elements,&
         num_outer_centers,&
         outer_center_coordinates,&
         outer_center_elements,&
         num_shells,&
         shell_centers,&
         shell_l_quantum_numbers,&
         shell_num_primitives,&
         primitive_exponents)

    grid => numgrid_get_grid(context)
    
!----------------------------------------------------------------------
! Number of grid points
!----------------------------------------------------------------------
    num_points=numgrid_get_num_points(context)

    return
    
  end subroutine initialise_intgrid

!######################################################################

  subroutine mo_cap_matrix(gam)

    use channels
    use constants
    use parameters
    use import_gamess
    use gamess_internal
    
    implicit none

    integer                               :: nao,i,j,k
    integer                               :: iatm,jatm,ish,jsh,ip,jp,&
                                             il,jl,bra,ket,icomp,jcomp,&
                                             inx,iny,inz,ipos,&
                                             jnx,jny,jnz,jpos
    real(dp)                              :: x,y,z,w,ix,iy,iz,jx,jy,jz,&
                                             aoi,aoj,ic,jc
    real(dp), dimension(:,:), allocatable :: sao,sao_grid
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! TEST: AO overlaps
!----------------------------------------------------------------------
    nao=gam%nbasis

    ! Analytic AO overlap matrix
    allocate(sao(nao,nao))
    call gamess_1e_integrals('AO OVERLAP',sao,gam,gam)

    ! Numerical AO overlap matrix
    !
    allocate(sao_grid(nao,nao))
    sao_grid=0.0d0
    
    ! Loop over bra AOs
    bra=0
    do iatm=1,gam%natoms               ! Loop over atoms
       do ish=1,gam%atoms(iatm)%nshell ! Loop over shells
          il=gam%atoms(iatm)%sh_l(ish) ! Angular momentum quantum no., l
          do icomp=1,gam_orbcnt(il)    ! Loop over components
                                       ! ({x,y,z},{xx,yy,...},etc) for
                                       ! the current l value

             bra=bra+1                 ! Bra counter


             ipos=ang_loc(il)+icomp-1  ! Position of the current AO in the ang_n*
                                       ! arrays

             inx=ang_nx(ipos)          ! nx
             iny=ang_ny(ipos)          ! ny
             inz=ang_nz(ipos)          ! nz

             ic=ang_c(ipos)
             
             ! Loop over ket AOs
             ket=0
             do jatm=1,gam%natoms               ! Loop over atoms
                do jsh=1,gam%atoms(jatm)%nshell ! Loop over shells
                   jl=gam%atoms(jatm)%sh_l(jsh) ! Angular momentum quantum no., l
                   do jcomp=1,gam_orbcnt(jl)    ! Loop over components
                                                ! ({x,y,z},{xx,yy,...},etc) for
                                                ! the current l value

                      ket=ket+1                 ! Ket counter


                      jpos=ang_loc(jl)+jcomp-1  ! Position of the current AO in the ang_n*
                                                ! arrays

                      jnx=ang_nx(jpos)          ! nx
                      jny=ang_ny(jpos)          ! ny
                      jnz=ang_nz(jpos)          ! nz

                      jc=ang_c(jpos)
                      
                      ! Calculation of the current integral
                      do k=1,num_points
                         
                         ! Coordinates and weight for the current
                         ! grid point
                         x=grid(k*4-3)
                         y=grid(k*4-2)
                         z=grid(k*4-1)
                         w=grid(k*4)
                         
                         ! Coordinate values relative to the atomic
                         ! centres (in Bohr)
                         ix=x-gam%atoms(iatm)%xyz(1)*ang2bohr
                         iy=y-gam%atoms(iatm)%xyz(2)*ang2bohr
                         iz=z-gam%atoms(iatm)%xyz(3)*ang2bohr
                         jx=x-gam%atoms(jatm)%xyz(1)*ang2bohr
                         jy=y-gam%atoms(jatm)%xyz(2)*ang2bohr
                         jz=z-gam%atoms(jatm)%xyz(3)*ang2bohr
                         
                         ! AO values
                         aoi=aoval(ix,iy,iz,inx,iny,inz,ic,iatm,ish,gam)
                         aoj=aoval(jx,jy,jz,jnx,jny,jnz,jc,jatm,jsh,gam)

                         ! Contribution to the integral
                         sao_grid(bra,ket)=sao_grid(bra,ket)+w*aoi*aoj
                         
                      enddo

                   enddo
                enddo
             enddo
             
          enddo
       enddo
    enddo
    

!    ! CHECK AO overlaps
!    print*,
!    do i=1,nao
!       do j=i,nao
!          print*,i,j,sao(i,j),sao_grid(i,j)
!       enddo
!    enddo
!    print*,
!    
!    STOP
    
    return
    
  end subroutine mo_cap_matrix

!######################################################################

  function aoval(x,y,z,nx,ny,nz,angc,iatom,ishell,gam)

    use constants
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
            * x**nx &
            * y**ny &
            * z**nz &
            * exp(-zeta*rsq)
        
    enddo

    return
    
  end function aoval
    
!######################################################################
  
  subroutine finalise_intgrid

    implicit none
    
    deallocate(center_coordinates)
    deallocate(center_elements)
    deallocate(shell_centers)
    deallocate(shell_l_quantum_numbers)
    deallocate(shell_num_primitives)
    deallocate(primitive_exponents)
    
    call numgrid_free_context(context)
    
  end subroutine finalise_intgrid
    
!######################################################################
  
end module capmod
