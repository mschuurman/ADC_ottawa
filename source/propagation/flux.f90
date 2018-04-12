!######################################################################
! fluxmod: routines for the calculation of flux expectation values
!######################################################################

module fluxmod

  implicit none

contains

!######################################################################
! adc1_flux_cap: calculation of the flux expectation value for a
!                TD-ADC(1) calculation using a CAP operator
!                analysis
!######################################################################
  
  subroutine adc1_flux_cap(matdim,psi,dtpsi,flux)

    use constants
    use parameters
    use iomod
    
    implicit none

    integer                               :: matdim
    integer                               :: unit,i,k
    real(d)                               :: flux,val1,val2,ener
    real(d), allocatable                  :: rvec(:)
    complex(d), dimension(matdim)         :: psi,dtpsi
    complex(d), dimension(:), allocatable :: oppsi
    complex(d), allocatable               :: vtmp1(:)
    complex(d), allocatable               :: vtmp2(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(oppsi(matdim))
    oppsi=czero

    allocate(vtmp1(matdim))
    allocate(vtmp2(matdim))
    vtmp1=czero
    vtmp2=czero
    
    if ((lprojcap.and.statenumber.gt.0) &
         .or.(statenumber.eq.0.and.iprojcap.eq.2)) then
       allocate(rvec(matdim))
       rvec=0.0d0
    endif
    
!----------------------------------------------------------------------
! (I) 2 Re < d Psi/dt | Theta | Psi >
!----------------------------------------------------------------------    
! Three pieces: (a) IS representation of the projector, Theta_IJ.
!               (b) Ground state projector matrix element, Theta_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, Theta_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    oppsi=czero

    ! Make a copy of the input vector to project against selected
    ! bound states
    vtmp1=psi

    ! First projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Open the ADC(1)/CIS vector file
       call freeunit(unit)
       open(unit,file='SCRATCH/initvecs',status='unknown',&
            access='sequential',form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,matdim-1
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          vtmp1(1:matdim-1)=vtmp1(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),psi(1:matdim-1))
       enddo
       ! Close the ADC(1)/CIS vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) vtmp1(matdim)=czero
    
    ! (a) IS-IS block
    !
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +matmul(thetaij,vtmp1(1:matdim-1))

    ! (b) Ground state-ground state element
    !
    oppsi(matdim)=oppsi(matdim)+theta00*vtmp1(matdim)
    
    ! (c) Ground state-IS block
    !
    oppsi(matdim)=oppsi(matdim) &
         +dot_product(theta0j(1:matdim-1),vtmp1(1:matdim-1))
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +theta0j(1:matdim-1)*vtmp1(matdim)

    ! Second projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Copy of oppsi
       vtmp2=oppsi
       ! Open the ADC(1)/CIS vector file
       call freeunit(unit)
       open(unit,file='SCRATCH/initvecs',status='unknown',&
            access='sequential',form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,matdim-1
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          oppsi(1:matdim-1)=oppsi(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),vtmp2(1:matdim-1))
       enddo
       ! Close the ADC(1)/CIS vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) oppsi(matdim)=czero
    
    val1=2.0d0*real(dot_product(dtpsi,oppsi))

!----------------------------------------------------------------------
! (II) 2 < Psi | W | Psi >
!----------------------------------------------------------------------    
! Three pieces: (a) IS representation of the CAP operator, W_IJ.
!               (b) Ground state CAP matrix element, W_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, W_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    oppsi=czero

    ! Make a copy of the input vector to project against selected
    ! bound states
    vtmp1=psi

    ! First projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Open the ADC(1)/CIS vector file
       call freeunit(unit)
       open(unit,file='SCRATCH/initvecs',status='unknown',&
            access='sequential',form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,matdim-1
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          vtmp1(1:matdim-1)=vtmp1(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),psi(1:matdim-1))
       enddo
       ! Close the ADC(1)/CIS vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) vtmp1(matdim)=czero
    
    ! (a) IS-IS block
    !
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +matmul(wij,vtmp1(1:matdim-1))

    ! (b) Ground state-ground state element
    !
    oppsi(matdim)=oppsi(matdim)+w00*vtmp1(matdim)
    
    ! (c) Ground state-IS block
    !
    oppsi(matdim)=oppsi(matdim) &
         +dot_product(w0j(1:matdim-1),vtmp1(1:matdim-1))
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +w0j(1:matdim-1)*vtmp1(matdim)

    ! Second projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Copy of oppsi
       vtmp2=oppsi
       ! Open the ADC(1)/CIS vector file
       call freeunit(unit)
       open(unit,file='SCRATCH/initvecs',status='unknown',&
            access='sequential',form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,matdim-1
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          oppsi(1:matdim-1)=oppsi(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),vtmp2(1:matdim-1))
       enddo
       ! Close the ADC(1)/CIS vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) oppsi(matdim)=czero
    
    val2=2.0d0*real(dot_product(psi,oppsi))
    
!----------------------------------------------------------------------
! Flux expectation value
!----------------------------------------------------------------------
    flux=val1+val2
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(oppsi)
    deallocate(vtmp1)
    deallocate(vtmp2)
    if ((lprojcap.and.statenumber.gt.0) &
         .or.(statenumber.eq.0.and.iprojcap.eq.2)) &
         deallocate(rvec)
    
    return
    
  end subroutine adc1_flux_cap

!######################################################################
! adc2_flux_cap: calculation of the flux expectation value for a
!                TD-ADC(2) calculation using a CAP operator
!                analysis
!######################################################################
  
  subroutine adc2_flux_cap(matdim,psi,dtpsi,flux)

    use constants
    use parameters
    use iomod
    use tdsemod, only: opxvec_ext
    
    implicit none

    integer                               :: matdim
    integer                               :: i,k,unit
    real(d)                               :: flux,val1,val2,ener
    real(d), allocatable                  :: rvec(:)
    complex(d), dimension(matdim)         :: psi,dtpsi
    complex(d), dimension(:), allocatable :: oppsi
    complex(d), allocatable               :: vtmp1(:)
    complex(d), allocatable               :: vtmp2(:)
    character(len=70)                     :: filename
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(oppsi(matdim))
    oppsi=czero

    allocate(vtmp1(matdim))
    allocate(vtmp2(matdim))
    vtmp1=czero
    vtmp2=czero
    
    if ((lprojcap.and.statenumber.gt.0) &
         .or.(statenumber.eq.0.and.iprojcap.eq.2)) then
       allocate(rvec(matdim))
       rvec=0.0d0
    endif

!----------------------------------------------------------------------
! (I) 2 Re < d Psi/dt | Theta | Psi >
!----------------------------------------------------------------------    
! Three pieces: (a) IS representation of the projector, Theta_IJ.
!               (b) Ground state projector matrix element, Theta_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, Theta_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    oppsi=czero

    ! Make a copy of the input vector to project against selected
    ! bound states
    vtmp1=psi

    ! First projection against selected bound states
    !
    ! Excited state contribution
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Open the ADC(2) vector file
       call freeunit(unit)
       open(unit,file=davname,status='old',access='sequential',&
            form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,davstates
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          vtmp1(1:matdim-1)=vtmp1(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),psi(1:matdim-1))
       enddo
       ! Close the ADC(2) vector file
       close(unit)
    endif
    !
    ! Ground state contribution
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) vtmp1(matdim)=czero
    
    ! (a) IS-IS block
    !
    ! Calculate W*v1
    filename='SCRATCH/theta'
    call opxvec_ext(matdim-1,vtmp1(1:matdim-1),&
         oppsi(1:matdim-1),filename,nbuf_theta)
    
    oppsi(1:matdim-1)=oppsi(1:matdim-1)+theta00*oppsi(1:matdim-1)

    ! (b) Ground state-ground state element
    !
    oppsi(matdim)=oppsi(matdim)+theta00*vtmp1(matdim)
    
    ! (c) Ground state-IS block
    !
    oppsi(matdim)=oppsi(matdim) &
         +dot_product(theta0j(1:matdim-1),vtmp1(1:matdim-1))
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +theta0j(1:matdim-1)*vtmp1(matdim)

    ! Second projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Copy of oppsi
       vtmp2=oppsi
       ! Open the ADC(2) vector file
       call freeunit(unit)
       open(unit,file=davname,status='old',access='sequential',&
            form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,davstates
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          oppsi(1:matdim-1)=oppsi(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),vtmp2(1:matdim-1))
       enddo
       ! Close the ADC(2) vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) oppsi(matdim)=czero
    
    val1=2.0d0*real(dot_product(dtpsi,oppsi))    

!----------------------------------------------------------------------
! (II) 2 < Psi | W | Psi >
!----------------------------------------------------------------------    
! Three pieces: (a) IS representation of the CAP operator, W_IJ.
!               (b) Ground state CAP matrix element, W_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, W_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    oppsi=czero

    ! Make a copy of the input vector to project against selected
    ! bound states
    vtmp1=psi

    ! First projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Open the ADC(2) vector file
       call freeunit(unit)
       open(unit,file=davname,status='old',access='sequential',&
            form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,davstates
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          vtmp1(1:matdim-1)=vtmp1(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),psi(1:matdim-1))
       enddo
       ! Close the ADC(2) vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) vtmp1(matdim)=czero
    
    ! (a) IS-IS block
    !
    filename='SCRATCH/cap'
    call opxvec_ext(matdim-1,vtmp1(1:matdim-1),&
         oppsi(1:matdim-1),filename,nbuf_cap)
    
    oppsi(1:matdim-1)=oppsi(1:matdim-1)+w00*oppsi(1:matdim-1)

    ! (b) Ground state-ground state element
    oppsi(matdim)=oppsi(matdim)+w00*vtmp1(matdim)

    ! (c) Ground state-IS block
    oppsi(matdim)=oppsi(matdim) &
         +dot_product(w0j(1:matdim-1),vtmp1(1:matdim-1))
    oppsi(1:matdim-1)=oppsi(1:matdim-1) &
         +w0j(1:matdim-1)*vtmp1(matdim)

    ! Second projection against selected bound states
    !
    ! Excited state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.gt.0).or.iprojcap.eq.2) then
       ! Copy of oppsi
       vtmp2=oppsi
       ! Open the ADC(2) vector file
       call freeunit(unit)
       open(unit,file=davname,status='old',access='sequential',&
            form='unformatted')
       ! Project the input vector onto the space orthogonal to
       ! the selected states
       do i=1,davstates
          read(unit) k,ener,rvec(1:matdim-1)
          if (ener.gt.projlim) exit
          if (iprojcap.eq.1.and.i.gt.statenumber) exit
          if (projmask(i).eq.0) cycle
          oppsi(1:matdim-1)=oppsi(1:matdim-1) &
               -rvec(1:matdim-1) &
               *dot_product(rvec(1:matdim-1),vtmp2(1:matdim-1))
       enddo
       ! Close the ADC(2) vector file
       close(unit)
    endif
    !
    ! Ground state contribution to the projector
    if ((iprojcap.eq.1.and.statenumber.eq.0) &
         .or.iprojcap.eq.2) oppsi(matdim)=czero
    
    val2=2.0d0*real(dot_product(psi,oppsi))
    
!----------------------------------------------------------------------
! Flux expectation value
!----------------------------------------------------------------------
    flux=val1+val2
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(oppsi)
    deallocate(vtmp1)
    deallocate(vtmp2)
    if ((lprojcap.and.statenumber.gt.0) &
         .or.(statenumber.eq.0.and.iprojcap.eq.2)) &
         deallocate(rvec)    
    
    return
    
  end subroutine adc2_flux_cap
    
!######################################################################
  
end module fluxmod
