!######################################################################
! capmod: routines for the calculation of the MO representation of
!         the CAP operator.
!
!         Two types of CAP are supported:
!
!         (1) Sigmoidal CAPs of the type used by Lopata et al
!             (see equations 6&7 in JCP, 145, 094105 (2016))
!
!         (2) Monomial-type CAPs
!             (see equations 6&7 in JCP, 115, 6853 (2001))
!
!         For sigmoidal CAPS, numerical integration is performed
!         using either Becke's partitioning scheme
!         (JCP 88, 2547 (1988)) or Gauss-Legendre quadrature.
!         Makes use of the NUMGRID libraries of Radovan Bast if
!         Becke's partitioning scheme is used.
!
!         For monomial CAPs, all matrix elements can be evaluated
!         analytically using the equations derived by Santra and
!         Cederbaum in JCP, 115, 6853 (2001).
!######################################################################

module capmod

  use numgrid
  use, intrinsic :: iso_c_binding, only: c_ptr
  
  implicit none  
  
  save
  
  ! Annoyingly, the gamess_internal module contains a variable
  ! named 'd', so we will use 'dp' here instead
  integer, parameter     :: dp=selected_real_kind(8)

  ! Conversion factors
  real(dp), parameter    :: ang2bohr=1.889725989d0

  ! Basis dimensions
  integer                :: npbas
  integer                :: nao

  ! Geometric centre of the molecule
  real(dp), dimension(3) :: gcent
  
  ! NUMGRID arrays and variables
  integer                :: min_num_angular_points
  integer                :: max_num_angular_points
  integer                :: num_points
  integer                :: num_centers
  integer                :: num_outer_centers
  integer                :: num_shells
  integer                :: num_primitives
  integer, allocatable   :: center_elements(:)
  integer, allocatable   :: shell_centers(:)
  integer, allocatable   :: shell_l_quantum_numbers(:)
  integer, allocatable   :: shell_num_primitives(:)
  integer, allocatable   :: outer_center_elements(:)
  real(dp), allocatable  :: outer_center_coordinates(:)
  real(dp), allocatable  :: primitive_exponents(:)
  real(dp)               :: radial_precision
  real(dp), allocatable  :: center_coordinates(:)
  real(dp), pointer      :: grid(:)
  type(c_ptr)            :: context

  ! Gauss-Legendre quadrature arrays and variables
  integer                :: ngp
  real(dp), allocatable  :: xabsc(:)
  real(dp), allocatable  :: weig(:)
  real(dp)               :: xa,xb,ya,yb,za,zb
  
  ! CAP arrays
  real(dp), allocatable  :: cap(:)
  real(dp), allocatable  :: cap_ao(:,:)
  real(dp), allocatable  :: vdwr(:)
  real(dp), parameter    :: dscale=3.5d0

  ! Monomial CAP arrays
  integer, allocatable   :: kmu(:,:)
  real(dp), dimension(3) :: cstrt
  real(dp), allocatable  :: amunu(:,:)
  real(dp), allocatable  :: Smunu(:,:)
  real(dp), allocatable  :: Rmunu(:,:,:)
  real(dp), allocatable  :: Rmu(:,:)
  real(dp), allocatable  :: alpha(:)
  real(dp), allocatable  :: Nmu(:)
  real(dp), allocatable  :: cap_pbas(:,:)
  
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
! Calculate the geometric centre of the molecule - this is used in
! various different schemes to evalaute CAP matrix elements, and so
! should be precalculated
!----------------------------------------------------------------------
    call calc_geomcent(gam)
    
!----------------------------------------------------------------------
! Calculate the MO CAP matrix elements
!----------------------------------------------------------------------
    if (igrid.eq.0) then
       ! Analytic evaluation of the MO CAP matrix elements. Note
       ! that this is only possible for a monomial CAP
       call cap_mobas_monomial_ana(gam,cap_mo)
    else
       ! Numerical integration of the MO CAP matrix elements
       if (igrid.eq.1) then
          ! Becke's integration scheme
          call cap_mobas_becke(gam,cap_mo)
       else if (igrid.eq.2) then
          ! Gauss-Legendre quadrature
          call cap_mobas_gauss(gam,cap_mo)
       endif
    endif
       
!----------------------------------------------------------------------
! Output timings
!----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Time taken:',tw2-tw1," s"
    
    return
    
  end subroutine cap_mobas

!######################################################################

  subroutine calc_geomcent(gam)

    use channels
    use parameters
    use timingmod
    use import_gamess

    implicit none

    integer             :: i,n,natom
    type(gam_structure) :: gam

!----------------------------------------------------------------------
! Geometric centre of the molecule (in Bohr)
!----------------------------------------------------------------------
    gcent=0.0d0
    natom=gam%natoms
    do n=1,natom
       do i=1,3
          gcent(i)=gcent(i)+gam%atoms(n)%xyz(i)*ang2bohr/natom
       enddo
    enddo

    return
    
  end subroutine calc_geomcent
    
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
    ! Integration boundaries
    xa=gcent(1)-0.5d0*gridpar(1)*ang2bohr
    xb=gcent(1)+0.5d0*gridpar(1)*ang2bohr
    ya=gcent(2)-0.5d0*gridpar(2)*ang2bohr
    yb=gcent(2)+0.5d0*gridpar(2)*ang2bohr
    za=gcent(3)-0.5d0*gridpar(3)*ang2bohr
    zb=gcent(3)+0.5d0*gridpar(3)*ang2bohr
    
!----------------------------------------------------------------------
! Calculate the MO representation of the CAP operator
!----------------------------------------------------------------------
    call mo_cap_matrix(gam,cap_mo)
    
    return
    
  end subroutine cap_mobas_gauss

!######################################################################
! cap_mobas_monomial_ana: Analytic evaluation of the MO representation
!                         of a monomial-type CAP operator.
!                         For a definition of the CAP used, see Eqs
!                         6 &7 in JCP, 115, 6853 (2001).
!                         Also see Eqs 9-19 of this paper for the
!                         working equations used to evaluate the
!                         matrix elements.
!######################################################################
  
  subroutine cap_mobas_monomial_ana(gam,cap_mo)

    use channels
    use iomod
    use parameters
    use import_gamess
    
    implicit none

    integer                               :: i,j,unit
    real(dp), dimension(:,:), allocatable :: cap_mo
    type(gam_structure)                   :: gam

!----------------------------------------------------------------------
! Fill in the van der Waals radius array
!----------------------------------------------------------------------
    call get_vdwr(gam)
    
!----------------------------------------------------------------------
! Precalculation of terms that appear a lot in the working equations
!----------------------------------------------------------------------
    call monomial_precalc(gam)
    
!----------------------------------------------------------------------
! Set up the CAP box: in each Cartesian direction, we take the start
! of the CAP to correspond to the furthest atom plus dscale times
! its van der Waals radius
!----------------------------------------------------------------------
    call get_cap_box_monomial(gam)
    
!----------------------------------------------------------------------
! Calculate the primitive representation of the CAP
!----------------------------------------------------------------------
    call primitive_monomial_cap_matrix(gam)

!----------------------------------------------------------------------
! Calculate the AO representation of the CAP
!----------------------------------------------------------------------
    call ao_monomial_cap_matrix(gam)

!----------------------------------------------------------------------
! Calculate the MO representation of the CAP
!----------------------------------------------------------------------
    ! MO representation of the CAP operator
    allocate(cap_mo(nbas,nbas))
    cap_mo=0.0d0

    ! Similarity transform the AO CAP matrix to yield the MO CAP matrix
    cap_mo=matmul(transpose(ao2mo),matmul(cap_ao,ao2mo))

    ! Multiplication by the CAP strength parameter
    cap_mo=cap_mo*capstr

    ! Symmetrisation - not really necessary but lets make sure that
    ! the MO CAP matrix is strictly Hermitian
    do i=1,nbas-1
       do j=i+1,nbas
          cap_mo(i,j)=cap_mo(j,i)
       enddo
    enddo

!----------------------------------------------------------------------
! For checking purposes, output the MO CAP matrix elements to file
!----------------------------------------------------------------------
!    call freeunit(unit)
!    open(unit,file='mocap.dat',form='formatted',status='unknown')
!    do i=1,nbas
!       do j=i,nbas
!          write(unit,'(2(i3,2x),ES15.7)') i,j,cap_mo(i,j)
!       enddo
!    enddo
!    close(unit)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(amunu)
    deallocate(Smunu)
    deallocate(Rmunu)
    deallocate(Rmu)
    deallocate(alpha)
    deallocate(kmu)
    deallocate(Nmu)
    deallocate(vdwr)
    deallocate(cap_pbas)
    deallocate(cap_ao)
    deallocate(ao2mo)
    
    return
    
  end subroutine cap_mobas_monomial_ana
    
!######################################################################

  subroutine monomial_precalc(gam)

    use channels
    use parameters
    use iomod
    use import_gamess
    use gamess_internal
    
    implicit none

    integer                :: i,j,l,ipr,jpr,lpr,mu,nu,n,m,mpr,il,&
                              ilpr,ipos
    real(dp)               :: zmu,znu
    real(dp)               :: ftmp
    type(gam_structure)    :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! No. primitives
    call get_npbas(gam)

    allocate(amunu(npbas,npbas))
    amunu=0.0d0

    allocate(Smunu(npbas,npbas))
    Smunu=0.0d0

    allocate(Rmunu(npbas,npbas,3))
    Rmunu=0.0d0

    allocate(Rmu(npbas,3))
    Rmu=0.0d0

    allocate(alpha(npbas))
    alpha=0.0d0

    allocate(kmu(npbas,3))
    kmu=0

    allocate(Nmu(npbas))
    Nmu=0.0d0

!----------------------------------------------------------------------
! R_mu - note that we take the atomic centres relative to the
!        geometric centre of the molecule as this makes the
!        construction of the monomial-type CAP easier
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1
                Rmu(mu,:)=gam%atoms(i)%xyz(1:3)*ang2bohr
                Rmu(mu,:)=Rmu(mu,:)-gcent(i)
             enddo
          enddo
       enddo
    enddo
    
!----------------------------------------------------------------------
! a_mu,nu
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1

                nu=0
                do ipr=1,gam%natoms
                   do jpr=1,gam%atoms(ipr)%nshell
                      ilpr=gam%atoms(ipr)%sh_l(jpr)
                      do mpr=1,gam_orbcnt(ilpr)
                         do lpr=gam%atoms(ipr)%sh_p(jpr),gam%atoms(ipr)%sh_p(jpr+1)-1
                            nu=nu+1
                            
                            amunu(mu,nu)=gam%atoms(i)%p_zet(l)&
                                 +gam%atoms(ipr)%p_zet(lpr)

                         enddo
                      enddo
                   enddo
                enddo
                      
             enddo
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! S_mu,nu
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1

                nu=0
                do ipr=1,gam%natoms
                   do jpr=1,gam%atoms(ipr)%nshell
                      ilpr=gam%atoms(ipr)%sh_l(jpr)
                      do mpr=1,gam_orbcnt(ilpr)
                         do lpr=gam%atoms(ipr)%sh_p(jpr),gam%atoms(ipr)%sh_p(jpr+1)-1
                            nu=nu+1
                            
                            zmu=gam%atoms(i)%p_zet(l)
                            znu=gam%atoms(ipr)%p_zet(lpr)

                            ftmp=dot_product(Rmu(mu,:)-Rmu(nu,:),&
                                 Rmu(mu,:)-Rmu(nu,:))
                      
                            Smunu(mu,nu)=exp(-zmu*znu*ftmp/amunu(mu,nu))
                            
                         enddo
                      enddo
                   enddo
                enddo
                      
             enddo
          enddo
       enddo
    enddo
    
!----------------------------------------------------------------------
! alpha_mu (exponents)
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1
                alpha(mu)=gam%atoms(i)%p_zet(l)
             enddo
          enddo
       enddo
    enddo
    
!----------------------------------------------------------------------
! R_mu,nu
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1

                nu=0
                do ipr=1,gam%natoms
                   do jpr=1,gam%atoms(ipr)%nshell
                      ilpr=gam%atoms(ipr)%sh_l(jpr)
                      do mpr=1,gam_orbcnt(ilpr)
                         do lpr=gam%atoms(ipr)%sh_p(jpr),gam%atoms(ipr)%sh_p(jpr+1)-1
                            nu=nu+1
                            
                            zmu=gam%atoms(i)%p_zet(l)
                            znu=gam%atoms(ipr)%p_zet(lpr)                     

                            Rmunu(mu,nu,:)=(zmu*Rmu(mu,:)+znu*Rmu(nu,:))/amunu(mu,nu)
                            
                         enddo
                      enddo
                   enddo
                enddo
                      
             enddo
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! k_mu,i
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             ipos=ang_loc(il)+m-1
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1
                
                kmu(mu,1)=ang_nx(ipos)
                kmu(mu,2)=ang_ny(ipos)
                kmu(mu,3)=ang_nz(ipos)

             enddo
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! N_mu (primitive normalisation constants)
!----------------------------------------------------------------------
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             ipos=ang_loc(il)+m-1
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1
                
                Nmu(mu)=ang_c(ipos)
                
             enddo
          enddo
       enddo
    enddo

    return
    
  end subroutine monomial_precalc

!######################################################################
! get_npbas: determines the total number of primitive basis functions
!######################################################################
  
  subroutine get_npbas(gam)

    use channels
    use parameters
    use iomod
    use import_gamess
    use gamess_internal
    
    implicit none

    integer             :: iatm,ish,il
    type(gam_structure) :: gam

    npbas=0
    do iatm=1,gam%natoms

       do ish=1,gam%atoms(iatm)%nshell

          il=gam%atoms(iatm)%sh_l(ish)

          npbas=npbas+gam_orbcnt(il)*(gam%atoms(iatm)%sh_p(ish+1) &
               -gam%atoms(iatm)%sh_p(ish))
           
       enddo
    enddo
    
    return
    
  end subroutine get_npbas

!######################################################################
! get_cap_box_monomial: determines the start of the rectangular CAP
!                       box used in a monomial CAP. For each Cartesian
!                       direction, we take the start of the CAP to
!                       correspond to the furthest atom plus its van
!                       der Waals radius multiplied by dscale
!######################################################################
  
  subroutine get_cap_box_monomial(gam)

    use channels
    use parameters
    use import_gamess
    use gamess_internal
    
    implicit none

    integer               :: k,i,j,l,mu,il,m
    integer, dimension(3) :: indx
    real(dp)              :: ftmp
    type(gam_structure)   :: gam

!----------------------------------------------------------------------
! Calculate the CAP box parameters. These are the c_i's of Eqs 6&7 in
! JCP, 115, 6853 (2001)
!
! Note that R_mu already has the geometric centre coordinates
! subtracted, so it is a bit simpler to work with these coordinates
!----------------------------------------------------------------------
    do k=1,3
       ftmp=0.0d0
       mu=0
       do i=1,gam%natoms

          do j=1,gam%atoms(i)%nshell
             il=gam%atoms(i)%sh_l(j)

             do m=1,gam_orbcnt(il)
             
                do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                   mu=mu+1
                   
                   if (abs(Rmu(mu,k))+dscale*vdwr(i).gt.ftmp) then
                      
                      ftmp=abs(Rmu(mu,k))+dscale*vdwr(i)
                      
                      cstrt(k)=abs(Rmu(mu,k))+dscale*vdwr(i)
                      
                   endif
                
                enddo
             enddo
          enddo
       enddo
    enddo

    return
    
  end subroutine get_cap_box_monomial
  
!######################################################################
! primitive_monomial_cap_matrix: calculates the primitive basis
!                                representation of a monomial CAP
!                                operator
!######################################################################
  
  subroutine primitive_monomial_cap_matrix(gam)

    use channels
    use parameters
    use iomod
    use import_gamess
    use gamess_internal
    
    implicit none

    integer               :: i,j,mu,nu,n
    real(dp), allocatable :: Xi(:,:,:)
    real(dp), allocatable :: Theta(:,:,:)
    real(dp)              :: Wi
    type(gam_structure)   :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Xi(npbas,npbas,3))
    Xi=0.0d0

    allocate(Theta(npbas,npbas,3))
    Theta=0.0d0

    allocate(cap_pbas(npbas,npbas))
    cap_pbas=0.0d0
    
!----------------------------------------------------------------------
! Temporary hardwiring of the CAP order
!----------------------------------------------------------------------
    n=2
    
!----------------------------------------------------------------------
! Calculate the terms Xi_mu,nu,i
!----------------------------------------------------------------------
    do i=1,3
       do mu=1,npbas
          do nu=1,npbas
             Xi(mu,nu,i)=Xival(mu,nu,i,n)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the terms Theta_mu,nu,j of Eq. 12 in JCP, 115, 6853 (2001)
!----------------------------------------------------------------------
    do j=1,3
       do mu=1,npbas
          do nu=1,npbas
             Theta(mu,nu,j)=Thetaval(mu,nu,j)
          enddo
       enddo
    enddo
    
!----------------------------------------------------------------------
! Calculate the primitive basis representation of the CAP operator.
! See Eq. 10 in JCP, 115, 6853 (2001)
!----------------------------------------------------------------------
    cap_pbas=0.0d0

    do i=1,3

       do mu=1,npbas
          do nu=1,npbas
             
             Wi=Nmu(mu)*Nmu(nu)*Smunu(mu,nu)*Xi(mu,nu,i)
             do j=1,3
                if (j.ne.i) Wi=Wi*Theta(mu,nu,j)
             enddo
             
             cap_pbas(mu,nu)=cap_pbas(mu,nu)+Wi
             
          enddo
       enddo

    enddo
          
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Xi)
    deallocate(Theta)
    
    return

  end subroutine primitive_monomial_cap_matrix

!######################################################################
! Xival: evaluates the terms Xi_mu,nu,i using Eq. 14 in
!        JCP, 115, 6853 (2001)
!######################################################################
  
  function Xival(mu,nu,i,n)

    use channels
    use parameters
    use iomod
    use math
    
    implicit none
    
    integer  :: mu,nu,i,n
    integer  :: rho,sigma
    real(dp) :: Xival
    real(dp) :: bcmu,bcnu,f1,f2,Imunu
    
    Xival=0.0d0
    
    do rho=0,kmu(mu,i)
       do sigma=0,kmu(nu,i)

          ! Binomial coefficients
          bcmu=mathbinomial(kmu(mu,i),rho)
          bcnu=mathbinomial(kmu(nu,i),sigma)

          ! (R_mu,nu,i - R_mu,i)^(k_mu,i-rho)
          f1=(Rmunu(mu,nu,i)-Rmu(mu,i))**(kmu(mu,i)-rho)

          ! (R_mu,nu,i - R_nu,i)^(k_nu,i-sigma)
          f2=(Rmunu(mu,nu,i)-Rmu(nu,i))**(kmu(nu,i)-sigma)

          ! I_mu,nu,i
          Imunu=Ival(mu,nu,i,cstrt(i),n,rho+sigma)

          Xival=Xival+bcmu*bcnu*f1*f2*Imunu
          
       enddo
    enddo
    
    return
    
  end function Xival

!######################################################################
! Thetaval: evaluates the terms Theta_mu,nu,i using Eq. 14 in
!           JCP, 115, 6853 (2001) and noting that
!
!           Theta_mu,nu,i = Xi_mu,nu,j(c=0,n=0)
!
!######################################################################
  
  function Thetaval(mu,nu,i)

    use channels
    use parameters
    use iomod
    use math
    
    implicit none
    
    integer  :: mu,nu,i
    integer  :: rho,sigma
    real(dp) :: Thetaval
    real(dp) :: bcmu,bcnu,f1,f2,Imunu
    
    Thetaval=0.0d0
    
    do rho=0,kmu(mu,i)
       do sigma=0,kmu(nu,i)

          ! Binomial coefficients
          bcmu=mathbinomial(kmu(mu,i),rho)
          bcnu=mathbinomial(kmu(nu,i),sigma)

          ! (R_mu,nu,i - R_mu,i)^(k_mu,i-rho)
          f1=(Rmunu(mu,nu,i)-Rmu(mu,i))**(kmu(mu,i)-rho)

          ! (R_mu,nu,i - R_nu,i)^(k_nu,i-sigma)
          f2=(Rmunu(mu,nu,i)-Rmu(nu,i))**(kmu(nu,i)-sigma)

          ! I_mu,nu,i
          Imunu=Ival(mu,nu,i,0.0d0,0,rho+sigma)
          
          Thetaval=Thetaval+bcmu*bcnu*f1*f2*Imunu
          
       enddo
    enddo
    
    return
    
  end function Thetaval
  
!######################################################################
! Ival: evaluates the terms I_mu,nu,i using Eq. 18 in
!       JCP, 115, 6853 (2001)
!######################################################################
  
  function Ival(mu,nu,i,c,n,kappa)

    use channels
    use parameters
    use iomod
    use math
    
    implicit none

    integer  :: mu,nu,i,n,kappa,tau
    real(dp) :: c,Ival
    real(dp) :: pre,f1,f2,Delta1,Delta2

    Ival=0.0d0

    do tau=0,n

       ! (-1)^(n-tau) * ('n choose tau') * a_mu,nu^-(kappa+tau+1)/2
       pre=(-1.0d0)**(n-tau) &
            * mathbinomial(n,tau) &
            * amunu(mu,nu)**(-0.5d0*(kappa+tau+1))

       ! (-1)^kappa * (c_i + R_mu,nu,i)^(n-tau)
       f1=(-1.0d0)**kappa &
            * (c+Rmunu(mu,nu,i))**(n-tau)

       ! (c_i - R_mu,nu,i)**(n-tau)
       f2=(c-Rmunu(mu,nu,i))**(n-tau)

       ! Delta(kappa+tau,a_mu,nu,c_i+R_mu,nu,i)
       Delta1=Deltaval(kappa+tau,&
                       amunu(mu,nu),&
                       c+Rmunu(mu,nu,i))
       
       ! Delta(kappa+tau,a_mu,nu,c_i-R_mu,nu,i)
       Delta2=Deltaval(kappa+tau,&
                       amunu(mu,nu),&
                       c-Rmunu(mu,nu,i))

       Ival=Ival+pre*(f1*Delta1 + f2*Delta2)
       
    enddo

    Ival=0.5d0*Ival
    
    return

  end function Ival

!######################################################################
! Deltaval: evaluates the Delta function of Eq. 19 in
!           JCP, 115, 6853 (2001)
!######################################################################
  
  function Deltaval(kappa,a,c)
    
    use channels
    use parameters
    use iomod
    use gammainc

    implicit none

    integer          :: kappa
    integer          :: ierr
    real(dp)         :: Deltaval,a,c
    real(dp)         :: g,ginc,xikappa
    character(len=1) :: ai

!----------------------------------------------------------------------
! Gamma function: Gamma((kappa+1)/2)
!----------------------------------------------------------------------
    g=gamma(dble(kappa+1)/2.0d0)

!----------------------------------------------------------------------
! Incomplete Gamma function: gamma((kappa+1)/2,ac^2)
!----------------------------------------------------------------------
    ! Incomplete gamma function
    ginc=g*gamma_inc_ratio(dble(kappa+1)/2.0d0,a*c**2)

!----------------------------------------------------------------------
! xi_kappa(c)
!----------------------------------------------------------------------
    if (c.ge.0.0d0) then
       xikappa=1.0d0
    else
       xikappa=(-1.0d0)**(kappa+1)
    endif

!----------------------------------------------------------------------
! Delta
!----------------------------------------------------------------------
    Deltaval=g-ginc*xikappa
    
    return
    
  end function Deltaval

!######################################################################
! ao_monomial_cap_matrix: calculates the AO representation of a
!                         monomial CAP operator from its primitive
!                         representation
!######################################################################
  
  subroutine ao_monomial_cap_matrix(gam)

    use channels
    use parameters
    use iomod
    use import_gamess
    use gamess_internal
    
    implicit none

    integer               :: i,j,il,m,l
    integer               :: q,mu
    real(dp), allocatable :: pbas2ao(:,:)
    type(gam_structure)   :: gam

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! No. AOs
    nao=gam%nbasis

    ! AO representation of the CAP operator
    allocate(cap_ao(nao,nao))
    cap_ao=0.0d0

    ! Primitive-to-AO transformation matrix
    allocate(pbas2ao(npbas,nao))
    pbas2ao=0.0d0
    
!----------------------------------------------------------------------
! Construct the primitive-to-AO transformation matrix
!----------------------------------------------------------------------
    ! Primitive counter
    mu=0

    ! AO counter
    q=0
    
    ! Loop over AOs
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             q=q+1
             
             ! Loop over the primitives that contribute to the
             ! current AO
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1

                ! Fill in the primitive-to-AO transformation matrix
                pbas2ao(mu,q)=gam%atoms(i)%p_c(l)

             enddo
                
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the AO representation of the CAP operator
!----------------------------------------------------------------------
    cap_ao=matmul(transpose(pbas2ao),matmul(cap_pbas,pbas2ao))
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(pbas2ao)
    
    return
    
  end subroutine ao_monomial_cap_matrix
    
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
    radial_precision=1.0e-12_d
    min_num_angular_points=170
    max_num_angular_points=974

!----------------------------------------------------------------------
! Atomic coordinates (in Bohr). Note that the gam_structure derived
! type holds the coordinates in Angstrom.
!----------------------------------------------------------------------
    ! Number of atoms
    num_centers=gam%natoms

    ! Atom coordinates
    allocate(center_coordinates(0:num_centers*3-1))
    center_coordinates=0.0d0
    
    k=0
    do i=1,gam%natoms
       do j=1,3
          k=k+1
          center_coordinates(k-1)=gam%atoms(i)%xyz(j)*ang2bohr
       enddo
    enddo

!----------------------------------------------------------------------
! Atom types
!----------------------------------------------------------------------
    allocate(center_elements(0:num_centers-1))
    center_elements=0
    do i=1,gam%natoms
       center_elements(i-1)=int(gam%atoms(i)%znuc)
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
    allocate(shell_centers(0:num_shells-1))
    shell_centers=0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_centers(k-1)=i
       enddo
    enddo

    ! Angular momentum quantum numbers for each shell
    allocate(shell_l_quantum_numbers(0:num_shells-1))
    shell_l_quantum_numbers=0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_l_quantum_numbers(k-1)=gam%atoms(i)%sh_l(j)
       enddo
    enddo

    ! Number of primitives for each shell
    allocate(shell_num_primitives(0:num_shells-1))
    shell_num_primitives=0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          k=k+1
          shell_num_primitives(k-1)=gam%atoms(i)%sh_p(j+1)-gam%atoms(i)%sh_p(j)
       enddo
    enddo

    ! Total number of primitives
    num_primitives=sum(shell_num_primitives)

    ! Primitive exponents
    allocate(primitive_exponents(0:num_primitives-1))
    primitive_exponents=0.0d0
    k=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
             k=k+1
             primitive_exponents(k-1)=gam%atoms(i)%p_zet(l)
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

    case(2) ! Monomial CAP

       call precalc_cap_monomial(gam)
       
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
    cap=0.0d0

    ! Loop over the grid points
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
! precalc_cap_monomial: Calculation of a monomial-type CAP at the
!                       grid points. The form of the CAP is given in
!                       equations 6 and 7 in JCP, 115, 6853 (2001).
!                       The CAP box parameters c_i are correspond, in
!                       each direction, to the position of the
!                       furthest atom plus its van der Waals radius
!                       multiplied by dscale
!######################################################################
  
  subroutine precalc_cap_monomial(gam)

    use channels
    use parameters
    use import_gamess
    use gamess_internal
    
    implicit none

    integer                :: i,n,mu,j,il,m,l,ord
    real(dp), dimension(3) :: x
    type(gam_structure)    :: gam
    
!----------------------------------------------------------------------
! Calculate the CAP box parameters.
! Note that the subroutine get_cap_box_monomial assumes that the Rmu
! array has already been filled in, so this must first be done.
!----------------------------------------------------------------------
    ! No. primitives
    call get_npbas(gam)

    ! R_mu - note that we take the atomic centres relative to the
    ! geometric centre of the molecule as this makes the
    ! construction of the monomial-type CAP easier
    allocate(Rmu(npbas,3))
    Rmu=0.0d0

    ! R_mu
    mu=0
    do i=1,gam%natoms
       do j=1,gam%atoms(i)%nshell
          il=gam%atoms(i)%sh_l(j)
          do m=1,gam_orbcnt(il)
             do l=gam%atoms(i)%sh_p(j),gam%atoms(i)%sh_p(j+1)-1
                mu=mu+1
                Rmu(mu,:)=gam%atoms(i)%xyz(1:3)*ang2bohr
                Rmu(mu,:)=Rmu(mu,:)-gcent(i)
             enddo
          enddo
       enddo
    enddo

    ! Calculate the CAP box parameters
    call get_cap_box_monomial(gam)
    
    deallocate(Rmu)

!----------------------------------------------------------------------
! Temporary hardwiring of the CAP order
!----------------------------------------------------------------------
    ord=2
    
!----------------------------------------------------------------------
! Evaluate the sigmoidal CAP value at each grid point
!----------------------------------------------------------------------
    cap=0.0d0

    ! Loop over the grid points
    do n=1,num_points

       ! Cartesian coordinates of the current grid point
       x(1)=grid(n*4-3)
       x(2)=grid(n*4-2)
       x(3)=grid(n*4-1)

       ! Loop over the x, y, and directions
       do i=1,3

          ! Contribution do the CAP
          if (abs(x(i)-gcent(i)).gt.cstrt(i)) then
             cap(n)=cap(n)+capstr*(abs(x(i))-cstrt(i))**ord
          endif
          
       enddo

    enddo
    
    return
    
  end subroutine precalc_cap_monomial
  
!######################################################################
! mo_cap_matrix: Numerical calculation of the MO representation of the
!                CAP operator. Additionally, for checking purposes,
!                the numerical AO overlap matrix is computed.
!######################################################################
  
  subroutine mo_cap_matrix(gam,cap_mo)

    use iomod
    use channels
    use parameters
    use import_gamess
    use gamess_internal
    use omp_lib
    
    implicit none

    integer                                 :: i,j,k,unit
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

             
             write(ilog,*) 'bra:',bra

             
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
! For checking purposes, output the MO CAP matrix elements to file
!----------------------------------------------------------------------
    !call freeunit(unit)
    !open(unit,file='mocap.dat',form='formatted',status='unknown')
    !do i=1,nbas
    !   do j=i,nbas
    !      write(unit,'(2(i3,2x),ES15.7)') i,j,cap_mo(i,j)
    !   enddo
    !enddo
    !close(unit)
    
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
       !write(ilog,*) i,sao(i,i),sao_grid(i,i)
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

          ! For monomial CAPs, we work with atomic coordinates
          ! taken relative to the goeometric centre of the
          ! molecule
          if (icap.eq.2) then
             ix=ix-gcent(1)
             iy=iy-gcent(2)
             iz=iz-gcent(3)
             jx=jx-gcent(1)
             jy=jy-gcent(2)
             jz=jz-gcent(3)
          endif          
          
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
    xm=0.5d0*(xb+xa)
    xl=0.5d0*(xb-xa)
    ym=0.5d0*(yb+ya)
    yl=0.5d0*(yb-ya)
    zm=0.5d0*(zb+za)
    zl=0.5d0*(zb-za)

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

       ! For monomial CAPs, we work with atomic coordinates
       ! taken relative to the goeometric centre of the
       ! molecule
       if (icap.eq.2) then
          ix=ix-gcent(1)
          jx=jx-gcent(1)
       endif
       
       do j=1,ngp
          ytmp=ym+yl*xabsc(j)
          iy=ytmp-gam%atoms(iatm)%xyz(1)*ang2bohr
          jy=ytmp-gam%atoms(jatm)%xyz(1)*ang2bohr

          if (icap.eq.2) then
             iy=iy-gcent(2)
             jy=jy-gcent(2)
          endif
          
          do k=1,ngp
             ztmp=zm+zl*xabsc(k)
             iz=ztmp-gam%atoms(iatm)%xyz(1)*ang2bohr
             jz=ztmp-gam%atoms(jatm)%xyz(1)*ang2bohr

             if (icap.eq.2) then
                iz=iz-gcent(3)
                jz=jz-gcent(3)
             endif
             
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
