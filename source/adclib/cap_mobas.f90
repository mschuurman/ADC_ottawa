!######################################################################
! capmod: routines for the calculation of the MO representation of
!         the CAP operator. Employs numerical integration performed
!         using either Becke's partitioning scheme
!         (JCP 88, 2547 (1988)) or Gauss-Legendre quadrature.
!         Makes use of the NUMGRID libraries of Radovan Bast if
!         Becke's partitioning scheme is used.
!######################################################################

module capmod

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

  ! Gauss-Legendre quadrature arrays and variables
  integer               :: ngp
  real(dp), allocatable :: xabsc(:)
  real(dp), allocatable :: weig(:)
  real(dp)              :: xi,xf,yi,yf,zi,zf
  
  ! CAP arrays
  real(dp), allocatable :: cap(:)
  real(dp), allocatable :: cap_ao(:,:)
  real(dp), allocatable :: vdwr(:)
  real(dp), parameter   :: dscale=3.5d0
  
contains
  
!######################################################################

  subroutine cap_mobas(gam,cap_mo)

    use channels
    use parameters
    use timingmod
    use import_gamess
    
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
    write(ilog,'(72a)') ('-',k=1,72)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Calculate the MO CAP matrix elements
!----------------------------------------------------------------------
    if (igrid.eq.1) then
       ! Becke's integration scheme
       call cap_mobas_becke(gam,cap_mo)
    else if (igrid.eq.2) then
       ! Gauss-Legendre quadrature
       call cap_mobas_gauss(gam,cap_mo)
    endif

!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine cap_mobas

!######################################################################

  subroutine cap_mobas_becke(gam,cap_mo)

    use channels
    use parameters
    use import_gamess
    
    implicit none

    real(dp), dimension(:,:), allocatable :: cap_mo
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! Fill in the van der Waals radius array
!----------------------------------------------------------------------
    call get_vdwr(gam)
    
!----------------------------------------------------------------------
! Set up the integration grid
!----------------------------------------------------------------------
    call initialise_intgrid(gam)

!----------------------------------------------------------------------
! Precalculation of the CAP at the grid points
!----------------------------------------------------------------------
    call precalc_cap(gam)
    
!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    call mo_cap_matrix(gam,cap_mo)

!----------------------------------------------------------------------
! Deallocate arrays and destroy the grid context
!----------------------------------------------------------------------
    call finalise_intgrid
    
    return
    
  end subroutine cap_mobas_becke

!######################################################################

  subroutine cap_mobas_gauss(gam,cap_mo)

    use channels
    use parameters
    use import_gamess
    
    implicit none

    integer                               :: i,n,natom
    real(dp), dimension(:,:), allocatable :: cap_mo
    real(dp), dimension(3)                :: cent
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! Fill in the van der Waals radius array
!----------------------------------------------------------------------
    call get_vdwr(gam)

!----------------------------------------------------------------------
! Calculate the Gauss-Legendre quadrature points and weights
!----------------------------------------------------------------------
    ngp=gridpar(4)
    allocate(xabsc(ngp))
    allocate(weig(ngp))
    
    call gauleg(ngp,xabsc,weig)

!----------------------------------------------------------------------
! Integration boundaries
!----------------------------------------------------------------------
    ! Geometric centre of the molecule (in Bohr)
    cent=0.0d0
    natom=gam%natoms
    do n=1,natom
       do i=1,3
          cent(i)=cent(i)+gam%atoms(n)%xyz(i)*ang2bohr/natom
       enddo
    enddo

    ! Integration boundaries
    xi=cent(1)-0.5d0*gridpar(1)*ang2bohr
    xf=cent(1)+0.5d0*gridpar(1)*ang2bohr
    yi=cent(2)-0.5d0*gridpar(2)*ang2bohr
    yf=cent(2)+0.5d0*gridpar(2)*ang2bohr
    zi=cent(3)-0.5d0*gridpar(3)*ang2bohr
    zf=cent(3)+0.5d0*gridpar(3)*ang2bohr
    
!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    call mo_cap_matrix(gam,cap_mo)
    
    return
    
  end subroutine cap_mobas_gauss
    
!######################################################################

  subroutine get_vdwr(gam)

    use channels
    use parameters
    use iomod
    use import_gamess
    
    implicit none

    integer             :: natom,i
    character(len=20)   :: name
    type(gam_structure) :: gam
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------   
    natom=gam%natoms
    allocate(vdwr(natom))
    vdwr=0.0d0

!----------------------------------------------------------------------
! Fill in the van der Waals radius array (units of Bohr)
! Van der Waals radii taken from Mantina et al., JPCA, 113, 5806 (2009)
!----------------------------------------------------------------------
    ! Van der Waals radii in Angstrom
    do i=1,natom

       name=gam%atoms(i)%name
       
       if (name.eq.'H') then
          vdwr(i)=1.10d0

       else if (name.eq.'He') then
          vdwr(i)=1.40d0

       else if (name.eq.'Li') then
          vdwr(i)=1.81d0

       else if (name.eq.'Be') then
          vdwr(i)=1.53d0

       else if (name.eq.'B') then
          vdwr(i)=1.92d0

       else if (name.eq.'C') then
          vdwr(i)=1.70d0

       else if (name.eq.'N') then
          vdwr(i)=1.55d0

       else if (name.eq.'O') then
          vdwr(i)=1.52d0

       else if (name.eq.'F') then
          vdwr(i)=1.47d0
                    
       else if (name.eq.'Ne') then
          vdwr(i)=1.54d0

       else if (name.eq.'Na') then
          vdwr(i)=2.27d0

       else if (name.eq.'Mg') then
          vdwr(i)=1.73d0

       else if (name.eq.'Al') then
          vdwr(i)=1.84d0

       else if (name.eq.'Si') then
          vdwr(i)=2.10d0

       else if (name.eq.'P') then
          vdwr(i)=1.80d0

       else if (name.eq.'S') then
          vdwr(i)=1.80d0

       else if (name.eq.'Cl') then
          vdwr(i)=1.75d0

       else if (name.eq.'Ar') then
          vdwr(i)=1.88d0

       else if (name.eq.'K') then
          vdwr(i)=2.75d0

       else if (name.eq.'Ca') then
          vdwr(i)=2.31d0

       else if (name.eq.'Ga') then
          vdwr(i)=1.87d0

       else if (name.eq.'Ge') then
          vdwr(i)=2.11d0

       else if (name.eq.'As') then
          vdwr(i)=1.85d0

       else if (name.eq.'Se') then
          vdwr(i)=1.90d0

       else if (name.eq.'Br') then
          vdwr(i)=1.83d0

       else if (name.eq.'Kr') then
          vdwr(i)=2.02d0

       else if (name.eq.'Rb') then
          vdwr(i)=3.03d0

       else if (name.eq.'Sr') then
          vdwr(i)=2.49d0

       else if (name.eq.'In') then
          vdwr(i)=1.93d0

       else if (name.eq.'Sn') then
          vdwr(i)=2.17d0

       else if (name.eq.'Sb') then
          vdwr(i)=2.06d0

       else if (name.eq.'Te') then
          vdwr(i)=2.06d0

       else if (name.eq.'I') then
          vdwr(i)=1.98d0

       else if (name.eq.'Xe') then
          vdwr(i)=2.16d0

       else if (name.eq.'Cs') then
          vdwr(i)=3.43d0

       else if (name.eq.'Ba') then
          vdwr(i)=2.68d0

       else if (name.eq.'Tl') then
          vdwr(i)=1.96d0

       else if (name.eq.'Pb') then
          vdwr(i)=2.02d0

       else if (name.eq.'Bi') then
          vdwr(i)=2.07d0

       else if (name.eq.'Po') then
          vdwr(i)=1.97d0

       else if (name.eq.'At') then
          vdwr(i)=2.02d0

       else if (name.eq.'Rn') then
          vdwr(i)=2.20d0

       else if (name.eq.'Fr') then
          vdwr(i)=3.48d0

       else if (name.eq.'Ra') then
          vdwr(i)=2.83d0

       else
          errmsg='Currently CAPs are not supported for the element '&
               //trim(name)
          call error_control
          
       endif
       
    enddo

    ! Convert to Bohr
    vdwr=vdwr*ang2bohr

    return
    
  end subroutine get_vdwr
  
!######################################################################

  subroutine initialise_intgrid(gam)

    use channels
    use parameters
    use import_gamess
    
    implicit none

    integer             :: i,j,k,l
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Grid parameters
!----------------------------------------------------------------------
    ! Grid parameters
    radial_precision=1.0e-20_d
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
    center_coordinates=0.0d0
    
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
    center_elements=0
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
    shell_centers=0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_centers(k)=i
       enddo
    enddo

    ! Angular momentum quantum numbers for each shell
    allocate(shell_l_quantum_numbers(num_shells))
    shell_l_quantum_numbers=0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_l_quantum_numbers(k)=gam%atoms(i)%sh_l(j)
       enddo
    enddo

    ! Number of primitives for each shell
    allocate(shell_num_primitives(num_shells))
    shell_num_primitives=0
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
    primitive_exponents=0.0d0
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
    ! Allocation of arrays that we do not need to use: depending on
    ! the compiler and compilation flags, we can run into problems
    ! when numgrid_generate_grid is called if these arrays are not
    ! allocated
    allocate(outer_center_elements(0))
    allocate(outer_center_coordinates(0))

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

  subroutine precalc_cap(gam)

    use channels
    use parameters
    use import_gamess
    use gamess_internal

    implicit none

    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(cap(num_points))
    cap=0.0d0
    
!----------------------------------------------------------------------
! Calculate the CAP values at the grid points
!----------------------------------------------------------------------
    select case(icap)

    case(1) ! Sigmoidal CAP

       call precalc_cap_sigmoidal(gam)
       
    end select
    
    return
    
  end subroutine precalc_cap

!######################################################################
! precalc_cap_sigmoidal: Calculation of a sigmoidal-type CAP at the
!                        grid points. The form of the CAP is given in
!                        equations 6 and 7 in JCP, 145, 094105 (2016).
!                        We choose the starting positions, R0, of the
!                        individual atom-centred CAPs to be 3.5 times
!                        the van der Waals radius of the atom. This
!                        value can be changed by altering the
!                        parameter dscale.
!######################################################################
  
  subroutine precalc_cap_sigmoidal(gam)

    use channels
    use parameters
    use import_gamess
    use gamess_internal

    implicit none

    integer                             :: natom,i,n
    real(dp)                            :: x,y,z,xa,ya,za,r,r0,&
                                           capval,capa
    real(dp)                            :: pival
    real(dp), dimension(:), allocatable :: acoo
    type(gam_structure)                 :: gam

    pival=4.0d0*atan(1.0d0)
    
!----------------------------------------------------------------------
! Atomic coordinates (in Bohr)
!----------------------------------------------------------------------
    natom=gam%natoms

    allocate(acoo(3*natom))
    acoo=0.0d0

    do n=1,natom
       do i=1,3
          acoo(n*3-3+i)=gam%atoms(n)%xyz(i)*ang2bohr
       enddo
    enddo
    
!----------------------------------------------------------------------
! Evaluate the sigmoidal CAP value at each grid point
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,num_points

       ! Cartesian coordinates of the current grid point
       x=grid(i*4-3)
       y=grid(i*4-2)
       z=grid(i*4-1)

       ! Loop over atoms
       capval=1e+12
       do n=1,natom

          ! Atomic coordinates
          xa=acoo(n*3-2)
          ya=acoo(n*3-1)
          za=acoo(n*3)

          ! Distance from the current grid point to the current atom
          r=sqrt((x-xa)**2+(y-ya)**2+(z-za)**2)
          
          ! CAP starting position
          r0=vdwr(n)*dscale

          ! Contribution to the total CAP value
          if (r.le.r0) then
             capa=0.0d0
          else if (r.gt.r0.and.r.lt.r0+capwid) then
             capa=capstr*(sin(pival*(r-r0)/(2.0d0*capwid)))**2
          else
             capa=capstr
          endif

          if (capa.lt.capval) capval=capa
          
       enddo

       cap(i)=capval
       
    enddo

    return
    
  end subroutine precalc_cap_sigmoidal
    
!######################################################################
! mo_cap_matrix: Numerical calculation of the MO representation of the
!                CAP operator. Additionally, for checking purposes,
!                the numerical AO overlap matrix is computed.
!######################################################################
  
  subroutine mo_cap_matrix(gam,cap_mo)

    use channels
    use parameters
    use import_gamess
    use gamess_internal
    use omp_lib
    
    implicit none

    integer                                 :: nao,i,j,k
    integer                                 :: iatm,jatm,ish,jsh,ip,jp,&
                                               il,jl,bra,ket,icomp,jcomp,&
                                               inx,iny,inz,ipos,&
                                               jnx,jny,jnz,jpos
    real(dp)                                :: iangc,jangc,&
                                               maxdiff,avdiff,trace
    real(dp), dimension(:,:), allocatable   :: sao,sao_grid,cap_mo,&
                                               sao_diff
    type(gam_structure)                     :: gam
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! No. AOs
    nao=gam%nbasis

    ! AO-to-MO transformation matrix
    if (.not.allocated(ao2mo)) then
       allocate(ao2mo(nao,nbas))
    endif
    ao2mo=0.0d0
    
    ! AO representation of the CAP operator
    allocate(cap_ao(nao,nao))
    cap_ao=0.0d0    
    
    ! MO representation of the CAP operator
    allocate(cap_mo(nbas,nbas))
    cap_mo=0.0d0
    
    ! Analytic AO overlap matrix
    allocate(sao(nao,nao))
    sao=0.0d0

    ! Numerical AO overlap matrix
    allocate(sao_grid(nao,nao))
    sao_grid=0.0d0

    ! Difference between the analytic and numerical AO overlap
    ! matrices
    allocate(sao_diff(nao,nao))
    sao_diff=0.0d0
    
!----------------------------------------------------------------------
! Analytic AO overlaps
!----------------------------------------------------------------------
    call gamess_1e_integrals('AO OVERLAP',sao,gam,gam)
    
!----------------------------------------------------------------------
! Calculate the AO representation of the CAP operator.
!
! For checking purposes, we also calculate the AO overlaps numerically
! and compare them to their analytic values. This way we can get some
! information about the quality of the integration grid.
!----------------------------------------------------------------------
    ! Loop over bra AOs
    bra=0
    do iatm=1,gam%natoms               ! Loop over atoms

       do ish=1,gam%atoms(iatm)%nshell ! Loop over shells

          il=gam%atoms(iatm)%sh_l(ish) ! Angular momentum quantum no., l

          do icomp=1,gam_orbcnt(il)    ! Loop over components
                                       ! ({x,y,z},{xx,yy,...},etc) for
                                       ! the current l value

             bra=bra+1                 ! Bra counter

             ipos=ang_loc(il)+icomp-1  ! Position of the current AO in
                                       ! the ang_n* arrays

             inx=ang_nx(ipos)          ! nx
             iny=ang_ny(ipos)          ! ny
             inz=ang_nz(ipos)          ! nz

             iangc=ang_c(ipos)

             
!             write(ilog,*) 'bra:',bra

             
             ! Loop over ket AOs
             ket=0
             do jatm=1,gam%natoms               ! Loop over atoms

                do jsh=1,gam%atoms(jatm)%nshell ! Loop over shells

                   jl=gam%atoms(jatm)%sh_l(jsh) ! Angular momentum quantum no., l

                   do jcomp=1,gam_orbcnt(jl)    ! Loop over components
                                                ! ({x,y,z},{xx,yy,...},etc) for
                                                ! the current l value

                      ket=ket+1                 ! Ket counter

                      jpos=ang_loc(jl)+jcomp-1  ! Position of the current AO in
                                                ! the ang_n* arrays

                      jnx=ang_nx(jpos)          ! nx
                      jny=ang_ny(jpos)          ! ny
                      jnz=ang_nz(jpos)          ! nz

                      jangc=ang_c(jpos)

                      ! Calculation of the current integral
                      if (igrid.eq.1) then
                         ! Becke's integration scheme
                         call calc_1integral_becke(cap_ao(bra,ket),&
                              sao_grid(bra,ket),gam,iatm,jatm,&
                              inx,iny,inz,jnx,jny,jnz,iangc,jangc,&
                              ish,jsh)
                      else if (igrid.eq.2) then
                         ! Gauss-Legendre quadrature
                         call calc_1integral_gauss(cap_ao(bra,ket),&
                              sao_grid(bra,ket),gam,iatm,jatm,&
                              inx,iny,inz,jnx,jny,jnz,iangc,jangc,&
                              ish,jsh)
                      endif
                      
                   enddo
                enddo
             enddo
             
          enddo
       enddo
    enddo    

!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    ! Retrieve the AO-to-MO transformation matrix
    ao2mo=gam%vectors(1:nao,1:nbas)

    ! Similarity transform the AO CAP matrix to yield the MO CAP matrix
    cap_mo=matmul(transpose(ao2mo),matmul(cap_ao,ao2mo))

!----------------------------------------------------------------------
! For checking purposes, output some information about the difference
! between the numerical and analytic AO overlap matrices and the
! trace of the numverical AO overlap matrix
!----------------------------------------------------------------------
    sao_diff=abs(sao-sao_grid)

    ! Maximum difference
    maxdiff=maxval(sao_diff)

    ! Average difference
    avdiff=0.0d0
    do i=1,nao
       do j=1,nao
          avdiff=avdiff+sao_diff(i,j)
       enddo
    enddo
    avdiff=avdiff/(nao**2)

    ! Trace of the numerical AO overlap matrix
    trace=0.0d0
    do i=1,nao
       trace=trace+sao_grid(i,i)
       write(ilog,*) i,sao(i,i),sao_grid(i,i)
    enddo

    write(ilog,'(/,2x,a,E15.7)') 'Trace of the numerical AO &
         overlap matrix:',trace
    
    write(ilog,'(/,2x,a)') 'AO overlap differences &
         (analytical - numerical):'
    write(ilog,'(2x,a,2x,E15.7)') 'Maximum:',maxdiff
    write(ilog,'(2x,a,2x,E15.7)') 'Average:',avdiff
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sao)
    deallocate(sao_grid)
    deallocate(cap_ao)
    deallocate(ao2mo)
    deallocate(sao_diff)
    
    return
    
  end subroutine mo_cap_matrix

!######################################################################
! calc_1integral_becke: Calculation of a single pair of AO CAP and
!                       overlap matrix elements using Becke's 
!                       integration scheme.
!######################################################################
  
  subroutine calc_1integral_becke(cap_ao,sao_grid,gam,iatm,jatm,inx,&
       iny,inz,jnx,jny,jnz,iangc,jangc,ish,jsh)

    use parameters
    use import_gamess
    use gamess_internal
    use omp_lib
    
    implicit none

    integer                              :: i,k
    integer                              :: nthreads,tid,npt
    integer, dimension(:,:), allocatable :: irange
    integer                              :: iatm,jatm,inx,iny,inz,&
                                            jnx,jny,jnz,ish,jsh
    real(dp)                             :: cap_ao,sao_grid
    real(dp)                             :: x,y,z,w,ix,iy,iz,jx,jy,jz,&
                                            aoi,aoj,iangc,jangc
    real(dp), dimension(:), allocatable  :: cap_ao_1thread,&
                                            sao_grid_1thread
    type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! Grid partitioning
    allocate(irange(nthreads,2))
    irange=0

    ! Working arrays
    allocate(cap_ao_1thread(nthreads))
    cap_ao_1thread=0.0d0

    allocate(sao_grid_1thread(nthreads))
    sao_grid_1thread=0.0d0
    
!-----------------------------------------------------------------------
! Partitioning of the grid points: one chunk per thread
!-----------------------------------------------------------------------
    npt=int(floor(real(num_points)/real(nthreads)))

    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=num_points
    
!-----------------------------------------------------------------------
! Calculate the current integral values
!-----------------------------------------------------------------------
    !$omp parallel do &
    !$omp& private(i,k,tid,x,y,z,w,ix,iy,iz,jx,jy,jz,aoi,aoj) &
    !$omp& shared(irange,sao_grid_1thread,cap_ao_1thread,grid,&
    !$omp& gam,cap,inx,iny,inz,iangc,ish,jnx,jny,jnz,jangc,jsh,&
    !$omp& iatm,jatm)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()
       
       do k=irange(tid,1),irange(tid,2)
          
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
          aoi=aoval(ix,iy,iz,inx,iny,inz,iangc,iatm,ish,gam)
          aoj=aoval(jx,jy,jz,jnx,jny,jnz,jangc,jatm,jsh,gam)
          
          ! Contribution to the AO overlap integral
          sao_grid_1thread(tid)=sao_grid_1thread(tid)+w*aoi*aoj
          
          ! Contribution to the AO CAP matrix element
          cap_ao_1thread(tid)=cap_ao_1thread(tid)+w*aoi*aoj*cap(k)
          
       enddo
                      
    enddo
    !$omp end parallel do

    ! Sum up the contributions from each thread
    sao_grid=0.0d0
    cap_ao=0.0d0
    do i=1,nthreads
       sao_grid=sao_grid+sao_grid_1thread(i)
       cap_ao=cap_ao+cap_ao_1thread(i)
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(irange)
    deallocate(cap_ao_1thread)
    deallocate(sao_grid_1thread)
    
    return
    
  end subroutine calc_1integral_becke

!######################################################################
! calc_1integral_gauss: Calculation of a single pair of AO CAP and
!                       overlap matrix elements using Gauss-Legendre
!                       quadrature.
!######################################################################
  
  subroutine calc_1integral_gauss(cap_ao,sao_grid,gam,iatm,jatm,inx,&
       iny,inz,jnx,jny,jnz,iangc,jangc,ish,jsh)

    use parameters
    use import_gamess
    use gamess_internal
    use omp_lib
    
    implicit none

    integer                              :: i,j,k
    integer                              :: nthreads,tid,npt
    integer, dimension(:,:), allocatable :: irange
    integer                              :: iatm,jatm,inx,iny,inz,&
                                            jnx,jny,jnz,ish,jsh
    real(dp)                             :: cap_ao,sao_grid
    real(dp)                             :: x,y,z,w,ix,iy,iz,jx,jy,jz,&
                                            aoi,aoj,iangc,jangc
    real(dp), dimension(:), allocatable  :: cap_ao_1thread,&
                                            sao_grid_1thread
    real(dp)                             :: xm,xl,ym,yl,zm,zl,xtmp,&
                                            ytmp,ztmp,stmp,captmp
    type(gam_structure)                  :: gam

!-----------------------------------------------------------------------
! Determine the no. threads
!-----------------------------------------------------------------------
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Calculate the matrix elements
!-----------------------------------------------------------------------
    xm=0.5d0*(xf+xi)
    xl=0.5d0*(xf-xi)
    ym=0.5d0*(yf+yi)
    yl=0.5d0*(yf-yi)
    zm=0.5d0*(zf+zi)
    zl=0.5d0*(zf-zi)

    stmp=0.0d0
    captmp=0.0d0
    
    !$omp parallel do &
    !$omp& private(i,j,k,xtmp,ytmp,ztmp,ix,iy,iz,jx,jy,jz,aoi,aoj) &
    !$omp& shared(xabsc,weig,inx,iny,inz,iangc,ish,jnx,jny,jnz,jangc,&
    !$omp& jsh,iatm,jatm) &
    !$omp& reduction(+:captmp,stmp)
    do i=1,ngp
       xtmp=xm+xl*xabsc(i)
       ix=xtmp-gam%atoms(iatm)%xyz(1)*ang2bohr
       jx=xtmp-gam%atoms(jatm)%xyz(1)*ang2bohr
       
       do j=1,ngp
          ytmp=ym+yl*xabsc(j)
          iy=ytmp-gam%atoms(iatm)%xyz(1)*ang2bohr
          jy=ytmp-gam%atoms(jatm)%xyz(1)*ang2bohr
          
          do k=1,ngp
             ztmp=zm+zl*xabsc(k)
             iz=ztmp-gam%atoms(iatm)%xyz(1)*ang2bohr
             jz=ztmp-gam%atoms(jatm)%xyz(1)*ang2bohr

             aoi=aoval(ix,iy,iz,inx,iny,inz,iangc,iatm,ish,gam)
             aoj=aoval(jx,jy,jz,jnx,jny,jnz,jangc,jatm,jsh,gam)

             stmp=stmp+weig(i)*weig(j)*weig(k)*aoi*aoj*xl*yl*zl
             captmp=captmp+weig(i)*weig(j)*weig(k)*capvalue(gam,xtmp,ytmp,ztmp)*aoi*aoj*xl*yl*zl
             
          enddo

       enddo

    enddo
    !$omp end parallel do

    sao_grid=stmp
    cap_ao=captmp
    
    return
    
  end subroutine calc_1integral_gauss
    
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
       aoval=aoval+coeff*angc*(x**nx)*(y**ny)*(z**nz)*exp(-zeta*rsq)
        
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
    deallocate(cap)
    deallocate(outer_center_elements)
    deallocate(outer_center_coordinates)
    
    call numgrid_free_context(context)

    return
    
  end subroutine finalise_intgrid

!######################################################################
! gauleg: Calculation of Gauss-Legendre quadrature points and weights.
!
! ngp: number of quadrature points/weights
! xabsc: quadrature points
! wei: quadrature weights
!
! Adapted from the gaussm3 code
!######################################################################
  
  subroutine gauleg(ngp,xabsc,weig)

    use constants
    
    implicit none

    integer                               :: i,j,m
    integer, intent(in)                   :: ngp
    real(dp)                              :: p1,p2,p3,pp,z,z1
    real(dp), dimension(ngp), intent(out) :: xabsc,weig
    real(dp), parameter                   :: eps=3.0d-15
    
!----------------------------------------------------------------------
! Number of desired roots.
! Note that the roots are symmetric and so we only need to compute
! half of them.
!----------------------------------------------------------------------
    m=(ngp+1)/2

!----------------------------------------------------------------------
! Compute the roots
!----------------------------------------------------------------------
    ! Loop over roots
    do i=1,m

       ! Approximation to the ith root
       z=cos(pi*(i-0.25d0)/(ngp+0.5d0))

       ! Refine the approximation to the ith root using Newton's
       ! method
100    p1=1.0d0
       p2=0.0d0

       ! Loop up the recurrence relation to get the Legendre
       ! polynomial evaluated at z
       do j=1,ngp
          p3=p2
          p2=p1
          p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
       enddo

       ! p1 is now the desired Legendre polynomial. We next compute pp,
       ! its derivative, by a standard relation involving also p2, the
       ! polynomial of one lower order
       pp=ngp*(z*p1-p2)/(z*z-1.0d0)
       z1=z
       z=z1-p1/pp
       if (dabs(z-z1).gt.eps) goto 100

       ! The roots are in the interval [-1,1]
       xabsc(i)=-z
       
       ! The roots are symmetric about the origin
       xabsc(ngp+1-i)=+z

       ! Calculate the ith weight
       weig(i)=2.0d0/((1.0d0-z*z)*pp*pp)

       ! Symmetric counterpart to the ith weight
       weig(ngp+1-i)=weig(i)
       
    enddo
    
    return
    
  end subroutine gauleg

!######################################################################  

  function capvalue(gam,x,y,z)

    use parameters
    use gamess_internal
    
    implicit none

    integer                             :: natom,i,n
    real(dp)                            :: x,y,z,capvalue
    real(dp)                            :: xa,ya,za,r,r0,currval,capa
    real(dp)                            :: pival
    real(dp), dimension(:), allocatable :: acoo
    type(gam_structure)                 :: gam
    
    pival=4.0d0*atan(1.0d0)

!----------------------------------------------------------------------
! Atomic coordinates (in Bohr)
!----------------------------------------------------------------------
    natom=gam%natoms

    allocate(acoo(3*natom))
    acoo=0.0d0

    do n=1,natom
       do i=1,3
          acoo(n*3-3+i)=gam%atoms(n)%xyz(i)*ang2bohr
       enddo
    enddo
    
!----------------------------------------------------------------------
! Evaluate the sigmoidal CAP value at the current point
!----------------------------------------------------------------------
    ! Initialisation
    currval=1e+12

    ! Loop over atoms
    do n=1,natom
       
       ! Atomic coordinates
       xa=acoo(n*3-2)
       ya=acoo(n*3-1)
       za=acoo(n*3)

       ! Distance from the current grid point to the current atom
       r=sqrt((x-xa)**2+(y-ya)**2+(z-za)**2)
       
       ! CAP starting position
       r0=vdwr(n)*dscale
       
       ! Contribution to the total CAP value
       if (r.le.r0) then
          capa=0.0d0
       else if (r.gt.r0.and.r.lt.r0+capwid) then
          capa=capstr*(sin(pival*(r-r0)/(2.0d0*capwid)))**2
       else
          capa=capstr
       endif
       
       if (capa.lt.currval) currval=capa
          
    enddo

    capvalue=currval

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------    
    deallocate(acoo)
    
    return
    
  end function capvalue
    
!######################################################################  
  
end module capmod
