!######################################################################
! Analytic evaluation of the MO representation of a monomial-type
! CAP operator. For a definition of the CAP used, see Eqs 6&7 in
! JCP, 115, 6853 (2001). Also see Eqs 9-19 of this paper for the
! working equations used to evaluate the matrix elements.
!######################################################################

module monomial_analytic

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

  subroutine monomial_ana(gam,cap_mo)

    use channels
    use parameters
    use misc, only: get_vdwr
    use import_gamess
    
    implicit none

    integer                        :: natom,i
    real(dp), dimension(nbas,nbas) :: cap_mo
    type(gam_structure)            :: gam

!----------------------------------------------------------------------
! Calculate the geometric centre of the molecule
!----------------------------------------------------------------------
    call calc_geomcent(gam)

!----------------------------------------------------------------------
! Fill in the van der Waals radius array
!----------------------------------------------------------------------
    natom=gam%natoms
    allocate(vdwr(natom))
    vdwr=0.0d0

    call get_vdwr(gam,vdwr,natom)

!----------------------------------------------------------------------
! Precalculation of terms that appear a lot in the working equations
!----------------------------------------------------------------------
    call monomial_precalc(gam)

!----------------------------------------------------------------------
! Set up the CAP box: if the CAP box has not been specified by the
! user, then in each Cartesian direction, we take the start of the CAP
! to correspond to the furthest atom plus dscale times its van der
! Waals radius
!----------------------------------------------------------------------
    if (boxpar(1).eq.0.0d0) then
       call get_cap_box_monomial(gam)
    else
       cstrt=boxpar
    endif

    write(ilog,'(2x,a,3(2x,F6.3))') 'CAP box parameters:',&
         (cstrt(i),i=1,3)

!----------------------------------------------------------------------
! Calculate the primitive representation of the CAP
!----------------------------------------------------------------------
    call primitive_monomial_cap_matrix(gam)

!----------------------------------------------------------------------
! Calculate the AO representation of the CAP
!----------------------------------------------------------------------
    call ao_cap_matrix(gam)

!----------------------------------------------------------------------
! Calculate the MO representation of the CAP
!----------------------------------------------------------------------
    call mo_cap_matrix(cap_mo)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    call finalise
    
    return
    
  end subroutine monomial_ana
    
!######################################################################
! calc_geomcent: calculation of the geometric centre of the molecule
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
! monomial_precalc: pre-calculation of a multitude of terms that enter
!                   into the analytic expressions for the primitive
!                   representation of a monomial CAP
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
! Calculate the terms Xi_mu,nu,i
!----------------------------------------------------------------------
    do i=1,3
       do mu=1,npbas
          do nu=1,npbas
             Xi(mu,nu,i)=Xival(mu,nu,i,capord)
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
! ao_cap_matrix: calculates the AO representation of the monomial CAP
!                operator from its primitive representation
!######################################################################

  subroutine ao_cap_matrix(gam)

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
    
  end subroutine ao_cap_matrix

!######################################################################
! mo_cap_matrix: calculates the MO representation of the monomial CAP
!                operator from its AO representation
!######################################################################
  
  subroutine mo_cap_matrix(cap_mo)
    
    use channels
    use iomod
    use parameters
    
    implicit none

    integer                        :: i,j
    real(dp), dimension(nbas,nbas) :: cap_mo
    real(dp), parameter            :: thrsh=1e-12_dp
    
!----------------------------------------------------------------------
! Similarity transform the AO CAP matrix to yield the MO CAP matrix
!----------------------------------------------------------------------
    cap_mo=matmul(transpose(ao2mo),matmul(cap_ao,ao2mo))

!----------------------------------------------------------------------
! Multiplication by the CAP strength parameter
!----------------------------------------------------------------------
    cap_mo=cap_mo*capstr

!----------------------------------------------------------------------
! Symmetry check
!----------------------------------------------------------------------
    do i=1,nbas-1
       do j=i+1,nbas
          if (abs(cap_mo(i,j)-cap_mo(j,i)).gt.thrsh) then
             errmsg='Error: the MO CAP matrix is not symmetric.'
             call error_control
          endif
          cap_mo(i,j)=cap_mo(j,i)
       enddo
    enddo

    return
    
  end subroutine mo_cap_matrix
    
!######################################################################

  subroutine finalise

    implicit none

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
    
    return
    
  end subroutine finalise
  
!######################################################################
  
end module monomial_analytic
