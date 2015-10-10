  module adc2dysonmod

    use channels

    contains

!#######################################################################

      subroutine adc2_dyson()

        use constants
        use parameters
        use fspace
        use misc
        use guessvecs

        implicit none

        integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
        integer                              :: ndim,ndims,ndimsf,&
                                                ndimf,ndimd
        real(d)                              :: e0
        real(d), dimension(:,:), allocatable :: eigvecf
        real(d), dimension(:), allocatable   :: eigvalf

!-----------------------------------------------------------------------
! Calculate the MP2 ground state energy and D2 diagnostic (if requested)
!-----------------------------------------------------------------------
        call run_mp2(e0)

!-----------------------------------------------------------------------  
! Initial space diagonalisation
!-----------------------------------------------------------------------  
        if (statenumber.gt.0) then
           write(ilog,'(/,2x,a,/)') 'The calculation of Dyson orbitals &
                for ionization from an excited state is not yet &
                supported'
           STOP
        endif

!-----------------------------------------------------------------------
! Final space diagonalisation (ionized states)
!-----------------------------------------------------------------------
        call calc_final_states(kpqf,ndimf,ndimsf,eigvecf,eigvalf)

!-----------------------------------------------------------------------
! Calculation of the expansion coefficients for the Dyson orbitals
! in the MO basis
!-----------------------------------------------------------------------
        call calc_dyscoeff(kpqf,eigvecf,eigvalf,ndimf,ndimsf)

        return

      end subroutine adc2_dyson

!#######################################################################

      subroutine run_mp2(e0)
        
        use constants
        use parameters
!        use diagnostics        
        use adc_ph, only: mp2

        implicit none

        integer           :: i
        real(d)           :: e0,d2
        character(len=60) :: atmp

!-----------------------------------------------------------------------
! Calculation of the MP2 correlation energy
!-----------------------------------------------------------------------
        call MP2(E_MP2)

!-----------------------------------------------------------------------
! Calculation of the D2 diagnostic
!-----------------------------------------------------------------------
!        if (ld2) call mp2_d2(d2)

!-----------------------------------------------------------------------
! Output results
!-----------------------------------------------------------------------
        e0=Ehf+E_MP2
 
        write(ilog,'(/,2x,90a)') ('*',i=1,90)

        write(ilog,'(2x,1a)') '*'
        atmp='* HF energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),Ehf

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 correlation energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),E_MP2

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),e0

        if (ld2) then
           write(ilog,'(2x,1a)') '*'
           atmp='* D2 diagnostic:'
           write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),d2
        endif

        write(ilog,'(2x,1a)') '*'
        write(ilog,'(2x,90a,/)') ('*',i=1,90)

        return

      end subroutine run_mp2

!#######################################################################

      subroutine calc_final_states(kpqf,ndimf,ndimsf,eigvecf,eigvalf)

        use constants
        use parameters

        implicit none
        
        integer, dimension(:,:), allocatable :: kpqf
        integer                              :: ndimf,ndimsf
        real(d), dimension(:,:), allocatable :: eigvecf
        real(d), dimension(:), allocatable   :: eigvalf

!-----------------------------------------------------------------------
! Determine the final subspace
!-----------------------------------------------------------------------
        call get_subspace_final(kpqf,ndimf,ndimsf)

!-----------------------------------------------------------------------
! Construct and diagonalise the IP-ADC(2)-s or IP-ADC(2)-x Hamiltonian 
! matrix
!-----------------------------------------------------------------------
        call diag_ipadc_mat(kpqf,ndimf,ndimsf,eigvecf,eigvalf)

        return

      end subroutine calc_final_states

!#######################################################################

      subroutine get_subspace_final(kpqf,ndimf,ndimsf)

        use constants
        use parameters
        use select_fano

        implicit none
        
        integer, dimension(:,:), allocatable :: kpqf
        integer                              :: ndimf,ndimsf

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(kpqf(7,0:nBas**2*4*nOcc**2))

!-----------------------------------------------------------------------
! Determine the final subspace
!-----------------------------------------------------------------------
        kpqf(:,:)=-1

        if (lcvsfinal) then
           call select_atom_isf_fakeip_cvs(kpqf(:,:))
           call select_atom_df_fakeip_cvs(kpqf(:,:),-1)
        else
           call select_atom_isf_fakeip(kpqf(:,:))
           call select_atom_df_fakeip(kpqf(:,:),-1)
        endif

        ndimsf=kpqf(1,0)
        ndimf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)

!-----------------------------------------------------------------------
! Output the final subspace information
!-----------------------------------------------------------------------
        write(ilog,*) 'ADC(2) FINAL Space dim',ndimf
        write(ilog,*) 'dimension of various FINAL configuration spaces'
        write(ilog,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
        
        return

      end subroutine get_subspace_final

!#######################################################################

      subroutine diag_ipadc_mat(kpqf,ndimf,ndimsf,eigvecf,eigvalf)

        use parameters
        use fspace
        use iomod, only: freeunit

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf
        real(d), dimension(:,:), allocatable      :: eigvecf
        real(d), dimension(:), allocatable        :: eigvalf
        character(len=120)                        :: msg

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        allocate(eigvecf(ndimf,ndimf))
        allocate(eigvalf(ndimf))

!-----------------------------------------------------------------------
! Calculate and diagonalise the IP-ADC(2) Hamiltonian matrix
!-----------------------------------------------------------------------
        if (method.eq.2) then
           msg='Diagonalising the IP-ADC(2)-s Hamiltonian matrix'
        else if (method.eq.3) then
           msg='Diagonalising the IP-ADC(2)-x Hamiltonian matrix'
        endif
        write(ilog,'(/,2x,a)') trim(msg)

        if (method.eq.2) then
           ! IP-ADC(2)-s
           call get_fspace_adc2_direct(ndimf,kpqf,eigvecf,eigvalf)
        else if (method.eq.3) then
           ! IP-ADC(2)-x
           call get_fspace_adc2e_direct(ndimf,kpqf,eigvecf,eigvalf)
        endif

        return

      end subroutine diag_ipadc_mat

!#######################################################################

      subroutine calc_dyscoeff(kpqf,eigvecf,eigvalf,ndimf,ndimsf)

        use constants
        use parameters
        use dysonmod
        use iomod, only: freeunit
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,i,j,&
                                                     a,b,n,inorm,icoeff
        real(d), dimension(ndimf,ndimf)           :: eigvecf
        real(d), dimension(ndimf)                 :: eigvalf
        real(d), dimension(:,:), allocatable      :: rhogs2
        real(d), dimension(:), allocatable        :: dyscoeff
        real(d)                                   :: norm
        
        write(ilog,'(/,2x,a,/)') 'Calculating the Dyson orbital &
             expansion coefficients'
        
!-----------------------------------------------------------------------
! Dyson orbitals for ionization from an excited state are not yet
! supported, so die here if this has been requested
!-----------------------------------------------------------------------
        if (statenumber.gt.0) then
           write(ilog,'(/,2x,a,/)') 'The calculation of Dyson orbitals &
                for ionization from an excited state is not yet &
                supported'
           STOP
        endif

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        ! Dyson orbital expansion coefficients
        allocate(dyscoeff(nbas))

        ! 2nd-order correction to the ground state density matrix
        allocate(rhogs2(nbas,nbas))

!-----------------------------------------------------------------------
! Pre-calculation of the 2nd-order correction to the ground state
! density matrix
!
! Note that this matrix is Hermitian, and so we only need to calculate
! the lower triangle
!-----------------------------------------------------------------------
        ! Occupied-occupied block
        do i=1,nocc
           do j=i,nocc
              rhogs2(i,j)=rhogs2_oo(i,j)
              rhogs2(j,i)=rhogs2(i,j)
           enddo
        enddo

        ! Unoccupied-unoccupied block
        do a=nocc+1,nbas
           do b=a,nbas
              rhogs2(a,b)=rhogs2_uu(a,b)
              rhogs2(b,a)=rhogs2(a,b)
           enddo
        enddo

        ! Occupied-unoccupied block
        do i=1,nocc
           do a=nocc+1,nbas
              rhogs2(i,a)=rhogs2_ou(i,a)
              rhogs2(a,i)=rhogs2(i,a)
           enddo
        enddo

!-----------------------------------------------------------------------
! Open output files
!-----------------------------------------------------------------------
        ! Dyson orbital norms
        call freeunit(inorm)
        open(inorm,file='dyson_norms',form='formatted',&
             status='unknown')

        ! Dyson orbital expansion coefficients
        call freeunit(icoeff)
        open(icoeff,file='dyson_coeffs',form='unformatted',&
             status='unknown')

!-----------------------------------------------------------------------        
! Calculation of the Dyson orbital expansion coefficients for all
! IP-ADC(2) states
!-----------------------------------------------------------------------
        do n=1,ndimf
           dyscoeff=0.0d0
           ! Coefficients for the occupied orbitals
           call get_dyscoeff_occ(rhogs2,dyscoeff,kpqf,eigvecf(:,n),&
                ndimf,ndimsf)
           ! Coefficients for the unoccupied orbitals
           call get_dyscoeff_unocc(rhogs2,dyscoeff,kpqf,eigvecf(:,n),&
                ndimf,ndimsf)           
           ! Output the Dyson orbital norm and coefficients
           norm=sqrt(dot_product(dyscoeff,dyscoeff))
           write(inorm,'(2(2x,F13.7))') eigvalf(n)*eh2ev,norm
        enddo

!-----------------------------------------------------------------------
! Close output files
!-----------------------------------------------------------------------
        close(inorm)
        close(icoeff)
        
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(rhogs2)
        
        return

      end subroutine calc_dyscoeff

!#######################################################################

      subroutine get_dyscoeff_occ(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)
        
        use constants        
        use parameters
        use vpqrsmod

        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,&
                                                     alpha,i,j,ilbl,&
                                                     jlbl,b,n,k,c
        real(d), dimension(nbas,nbas)             :: rhogs2
        real(d), dimension(nbas)                  :: dyscoeff
        real(d), dimension(ndimf)                 :: vec
        real(d), dimension(:,:), allocatable      :: chi,zeta
        real(d)                                   :: delta_ijab,ftmp
        
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
        ! Note that these arrays only need to be of
        ! dimension (nocc,nvirt), but for simplicity they are of
        ! dimension (nbas,nbas), as this requires only a negligible
        ! extra amount of memory 
        allocate(chi(nbas,nbas))
        allocate(zeta(nbas,nbas))
        
!-----------------------------------------------------------------------
! Index of the 'continuum' orbital
!-----------------------------------------------------------------------
        alpha=ifakeorb

!-----------------------------------------------------------------------
! Zeroth-order contribution
!-----------------------------------------------------------------------
        do i=1,ndimsf
           ilbl=kpqf(3,i)
           dyscoeff(ilbl)=vec(i)
        enddo

!-----------------------------------------------------------------------
! Second-oder contribution
! Note that the first-order contribution vanishes
!-----------------------------------------------------------------------
        ! Term C
        ! (i) Pre-calculate the chi and zeta terms
        chi=0.0d0
        zeta=0.0d0
        do j=1,nocc
           do b=nocc+1,nbas
              do n=1,ndimsf
                 k=kpqf(3,n)
                 c=kpqf(5,n)
                 chi(j,b)=chi(j,b)+vec(n)*(vpqrs(b,j,c,k)-vpqrs(c,j,b,k))
                 zeta(j,b)=zeta(j,b)+vec(n)*(vpqrs(c,j,b,k)-vpqrs(b,j,c,k))
              enddo
           enddo
        enddo        
        ! (ii) Contribution of the C terms
        do i=1,nocc
           do j=1,nocc
              do b=nocc+1,nbas
                 delta_ijab=1.0d0/(e(alpha)+e(b)-e(j)-e(j))
                 ftmp=vpqrs(i,alpha,j,b)*chi(j,b)+vpqrs(i,b,j,alpha)*zeta(j,b)
                 ftmp=delta_ijab*ftmp
                 dyscoeff(i)=dyscoeff(i)+ftmp
              enddo
           enddo
        enddo

        ! Term D
        do i=1,ndimsf
           ilbl=kpqf(3,i)
           dyscoeff(ilbl)=dyscoeff(ilbl)-0.5d0*rhogs2(alpha,alpha)*vec(i)
        enddo

        ! Term E
        do i=1,nocc
           do j=1,ndimsf
              jlbl=kpqf(3,j)
              dyscoeff(i)=dyscoeff(i)+0.5d0*rhogs2(i,jlbl)*vec(j)
           enddo
        enddo

!! CHECK COEFFICIENTS
!        ftmp=0.0d0
!        do i=1,nocc
!           ftmp=ftmp+dyscoeff(i)**2
!           print*,i,dyscoeff(i)
!        enddo
!        ftmp=sqrt(ftmp)
!        print*,"Norm (occ):",ftmp

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
        deallocate(chi)
        deallocate(zeta)
        
        return

      end subroutine get_dyscoeff_occ

!#######################################################################

      subroutine get_dyscoeff_unocc(rhogs2,dyscoeff,kpqf,vec,ndimf,&
           ndimsf)

        use constants
        use parameters
        use vpqrsmod
        
        implicit none

        integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
        integer                                   :: ndimf,ndimsf,&
                                                     alpha,b,i,ilbl,&
                                                     j,c,n
        real(d), dimension(nbas,nbas)             :: rhogs2
        real(d), dimension(nbas)                  :: dyscoeff
        real(d), dimension(ndimf)                 :: vec
        real(d)                                   :: delta_ijbc,ftmp
        
!-----------------------------------------------------------------------
! Index of the 'continuum' orbital
!-----------------------------------------------------------------------
        alpha=ifakeorb
        
!-----------------------------------------------------------------------
! Second-order contributions
! Note that the zeroth and first-order contributions vanish
!-----------------------------------------------------------------------
        ! Term F
        do b=nocc+1,nbas
           do i=1,nocc
              ilbl=kpqf(3,i)
              dyscoeff(b)=dyscoeff(b)+vec(i)*rhogs2(ilbl,b)
           enddo
        enddo

        ! Term G
        do b=nocc+1,nbas
           do n=ndimsf+1,ndimf
              i=kpqf(3,n)
              j=kpqf(4,n)
              if (kpqf(5,n).eq.alpha) then
                 c=kpqf(6,n)
              else
                 c=kpqf(5,n)
              endif
              delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))
              ftmp=vec(n)*(2.0d0*vpqrs(b,i,c,j)-vpqrs(c,i,b,j))
              ftmp=delta_ijbc*ftmp
              dyscoeff(b)=dyscoeff(b)+ftmp
           enddo
        enddo
        
        return
        
      end subroutine get_dyscoeff_unocc
        
!#######################################################################

  end module adc2dysonmod
