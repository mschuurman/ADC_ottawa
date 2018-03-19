!#######################################################################
! tsdemod: common routines used in electronic wavepacket propagation
!          calculations.
!#######################################################################

module tdsemod

  use constants
  
  save

  integer                            :: mxbl,nrecord,nmult
  integer, dimension(:), allocatable :: iindx,jindx
  real(d), dimension(:), allocatable :: hon,hoff
  logical                            :: hincore
  
contains

!#######################################################################
! load_hamiltonian: loads the non-zero elements of the Hamiltonian
!                   matrix into memory
!#######################################################################

  subroutine load_hamiltonian(fileon,fileoff,matdim,noffdiag)

    use iomod
    
    implicit none

    integer, intent(in)                :: matdim
    integer*8, intent(in)              :: noffdiag
    integer                            :: unit,count,k,nlim
    integer, dimension(:), allocatable :: itmp,jtmp
    real(d), dimension(:), allocatable :: ftmp
    character(len=*), intent(in)       :: fileon,fileoff

!-----------------------------------------------------------------------
! Get the next free file unit
!-----------------------------------------------------------------------
    call freeunit(unit)
    
!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    if (allocated(hon)) deallocate(hon)
    allocate(hon(matdim))

    open(unit,file=fileon,status='old',access='sequential',&
         form='unformatted')

    read(unit) mxbl,nrecord
    read(unit) hon
    
    close(unit)

!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    if (allocated(hoff)) deallocate(hoff)
    if (allocated(iindx)) deallocate(iindx)
    if (allocated(jindx)) deallocate(jindx)
    allocate(hoff(noffdiag))
    allocate(iindx(noffdiag))
    allocate(jindx(noffdiag))

    allocate(ftmp(mxbl))
    allocate(itmp(mxbl))
    allocate(jtmp(mxbl))

    open(unit,file=fileoff,status='old',access='sequential',&
           form='unformatted')

    count=0
    do k=1,nrecord
       read(unit) ftmp(:),itmp(:),jtmp(:),nlim
       hoff(count+1:count+nlim)=ftmp(1:nlim)
       iindx(count+1:count+nlim)=itmp(1:nlim)
       jindx(count+1:count+nlim)=jtmp(1:nlim)
       count=count+nlim
    enddo
      
    close(unit)
    
    deallocate(ftmp,itmp,jtmp)

    return
    
  end subroutine load_hamiltonian

!#######################################################################
! deallocate_hamiltonian: deallocates the hamiltonian arrays
!#######################################################################

  subroutine deallocate_hamiltonian

    implicit none
    
    if (allocated(hon)) deallocate(hon)
    if (allocated(hoff)) deallocate(hoff)
    if (allocated(iindx)) deallocate(iindx)
    if (allocated(jindx)) deallocate(jindx)
    
    return
    
  end subroutine deallocate_hamiltonian
  
!#######################################################################
! matxvec: calculates v2 = -H.v1 (Note the minus sign!)
!          (used in imaginary time wavepacket propagations)
!#######################################################################

  subroutine matxvec(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
  
    implicit none

    integer                    :: matdim
    integer*8                  :: noffdiag
    real(d), dimension(matdim) :: v1,v2

    if (hincore) then
       call matxvec_incore(matdim,noffdiag,v1,v2)
    else
       call matxvec_ext(matdim,noffdiag,v1,v2)
    endif

    ! Update nmult
    nmult=nmult+1
    
    return
    
  end subroutine matxvec

!#######################################################################
! matxvec_incore: calculates v2 = -H.v1 using a Hamiltonian stored
!                 in-core
!#######################################################################

  subroutine matxvec_incore(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
    
    implicit none

    integer                                :: matdim,maxbl,nrec,nlim,&
                                              ion,i,j,k,l,nthreads,tid
    integer*8                              :: noffdiag,npt
    integer*8, dimension(:,:), allocatable :: irange
    real(d), dimension(matdim)             :: v1,v2
    real(d), dimension(:,:), allocatable   :: tmpvec
    character(len=70)                      :: fileon,fileoff

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(irange(nthreads,2))
    allocate(tmpvec(matdim,nthreads))

!-----------------------------------------------------------------------
! Partitioning of the off-diagonal elements: one chunk per thread
!-----------------------------------------------------------------------
    npt=int(floor(real(noffdiag)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffdiag

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=-hon(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=0.0d0

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(iindx(k),tid)=tmpvec(iindx(k),tid)&
               -hoff(k)*v1(jindx(k))
          tmpvec(jindx(k),tid)=tmpvec(jindx(k),tid)&
               -hoff(k)*v1(iindx(k))
       enddo
          
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo
       
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(irange)
    deallocate(tmpvec)
    
    return
    
  end subroutine matxvec_incore
  
!#######################################################################
! matxvec_ext: calculates v2 = -H.v1 using an externally stored
!              Hamiltonian
!#######################################################################
  
  subroutine matxvec_ext(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
  
    implicit none

    integer                            :: matdim,maxbl,nrec,nlim,&
                                          ion,ioff,k,l
    integer*8                          :: noffdiag
    integer, dimension(:), allocatable :: indxi,indxj
    real(d), dimension(matdim)         :: v1,v2
    real(d), dimension(:), allocatable :: hii,hij
    character(len=70)                  :: fileon,fileoff

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
    ! Filenames
    if (hamflag.eq.'i') then
       fileon='SCRATCH/hmlt.diai'
       fileoff='SCRATCH/hmlt.offi'
    else if (hamflag.eq.'f') then
       fileon='SCRATCH/hmlt.diac'
       fileoff='SCRATCH/hmlt.offc'
    endif

    ! On-diagonal elements
    call freeunit(ion)
    open(ion,file=fileon,status='old',access='sequential',&
         form='unformatted')
    
    ! Off-diagonal elements
    call freeunit(ioff)
    open(ioff,file=fileoff,status='old',access='sequential',&
         form='unformatted')
  
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    rewind(ion)
    read(ion) maxbl,nrec

    allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
    allocate(hii(matdim))

!-----------------------------------------------------------------------
! On-diagonal elements
!-----------------------------------------------------------------------
    read(ion) hii
  
    do k=1,matdim
       v2(k)=-hii(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
    rewind(ioff)

    do k=1,nrec
       read(ioff) hij(:),indxi(:),indxj(:),nlim
       do l=1,nlim
          v2(indxi(l))=v2(indxi(l))-hij(l)*v1(indxj(l))
          v2(indxj(l))=v2(indxj(l))-hij(l)*v1(indxi(l))
       enddo
    enddo
  
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(hii,hij,indxi,indxj)
  
!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(ion)
    close(ioff)
    
    return
    
  end subroutine matxvec_ext

!#######################################################################
! matxvec_treal: calculates v2 = -iH.v1 
!                (used in real time wavepacket propagations)
!#######################################################################

  subroutine matxvec_treal(time,matdim,noffdiag,v1,v2)

    use constants

    implicit none
    
    integer                       :: matdim
    integer*8                     :: noffdiag
    real(d)                       :: time
    complex(d), dimension(matdim) :: v1,v2

    if (hincore) then
       call matxvec_treal_incore(matdim,noffdiag,v1,v2)
    else
       call matxvec_treal_ext(matdim,noffdiag,v1,v2)
    endif

    ! Update nmult
    nmult=nmult+1
    
    return

  end subroutine matxvec_treal

!#######################################################################
! matxvec_treal_incore: calculates v2 = -iH.v1 using a Hamiltonian
!                       stored in-core
!#######################################################################

  subroutine matxvec_treal_incore(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
  
    implicit none

    integer                                 :: matdim,maxbl,nrec,nlim,&
                                               ion,i,j,k,l,nthreads,tid
    integer*8                               :: noffdiag,npt
    integer*8, dimension(:,:), allocatable  :: irange
    complex(d), dimension(matdim)           :: v1,v2
    complex(d), dimension(:,:), allocatable :: tmpvec
    character(len=70)                       :: fileon,fileoff

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(irange(nthreads,2))
    allocate(tmpvec(matdim,nthreads))
    
!-----------------------------------------------------------------------
! Partitioning of the off-diagonal elements: one chunk per thread
!-----------------------------------------------------------------------
    npt=int(floor(real(noffdiag)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffdiag

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=-ci*hon(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=czero

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(iindx(k),tid)=tmpvec(iindx(k),tid)&
               -ci*hoff(k)*v1(jindx(k))
          tmpvec(jindx(k),tid)=tmpvec(jindx(k),tid)&
               -ci*hoff(k)*v1(iindx(k))
       enddo
          
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo
       
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(irange)
    deallocate(tmpvec)
    
    return

  end subroutine matxvec_treal_incore
    
!#######################################################################
! matxvec_treal_ext: calculates v2 = -iH.v1 using an externally
!                    stored Hamiltonian
!#######################################################################
  
  subroutine matxvec_treal_ext(matdim,noffdiag,v1,v2)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
  
    implicit none

    integer                                 :: matdim,maxbl,nrec,nlim,&
                                               ion,i,j,k,l,nthreads,tid
    integer*8                               :: noffdiag
    integer, dimension(:), allocatable      :: indxi,indxj,ioff
    real(d), dimension(:), allocatable      :: hii,hij
    complex(d), dimension(matdim)           :: v1,v2
    complex(d), dimension(:,:), allocatable :: tmpvec
    character(len=70)                       :: fileon,fileoff

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
    ! On-diagonal file
    if (hamflag.eq.'i') then
       fileon='SCRATCH/hmlt.diai'
    else
       fileon='SCRATCH/hmlt.diac'
    endif
    call freeunit(ion)
    open(ion,file=fileon,status='old',access='sequential',&
         form='unformatted')
    
    ! Off-diagonal files
    allocate(ioff(nthreads))
    do i=1,nthreads
       call freeunit(ioff(i))
       if (hamflag.eq.'i') then
          fileoff='SCRATCH/hmlt.offi.'
       else
          fileoff='SCRATCH/hmlt.offc.'
       endif
       k=len_trim(fileoff)+1
       if (i.lt.10) then
          write(fileoff(k:k),'(i1)') i
       else
          write(fileoff(k:k+1),'(i2)') i
       endif
       open(unit=ioff(i),file=fileoff,status='unknown',&
            access='sequential',form='unformatted')
    enddo

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    rewind(ion)
    read(ion) maxbl,nrec

    allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))
    allocate(hii(matdim))
    allocate(tmpvec(matdim,nthreads))

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements
!-----------------------------------------------------------------------
    read(ion) hii
  
    do k=1,matdim
       v2(k)=-ci*hii(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contributoion from the off-diagonal elements
!-----------------------------------------------------------------------
    tmpvec=czero

    !$omp parallel do &
    !$omp& private(i,k,l,tid,hij,indxi,indxj,nlim) &
    !$omp& shared(tmpvec,v1,nrec_omp)
    do i=1,nthreads
       tid=1+omp_get_thread_num()
       do k=1,nrec_omp(i)
          read(ioff(i)) hij(:),indxi(:),indxj(:),nlim
          do l=1,nlim
             tmpvec(indxi(l),tid)=tmpvec(indxi(l),tid) &
                  -ci*hij(l)*v1(indxj(l))
             tmpvec(indxj(l),tid)=tmpvec(indxj(l),tid) &
                  -ci*hij(l)*v1(indxi(l))
          enddo
       enddo
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(ion)

    do i=1,nthreads
       close(ioff(i))
    enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(hii,hij,indxi,indxj)
    deallocate(tmpvec)
    deallocate(ioff)
    
    return

  end subroutine matxvec_treal_ext

!#######################################################################

  subroutine matxvec_treal_laser(time,matdim,noffdiag,v1,v2)

    use constants
    use parameters
    
    implicit none
    
    integer                               :: matdim,i,j
    integer*8                             :: noffdiag
    real(d)                               :: time
    real(d), dimension(3)                 :: Et
    complex(d), dimension(matdim)         :: v1,v2
    complex(d), dimension(:), allocatable :: opxv1
    character(len=70)                     :: filename
    character(len=1), dimension(3)        :: acomp

    acomp=(/ 'x','y','z' /)
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    v2=czero
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(opxv1(matdim))
    opxv1=czero
    
!----------------------------------------------------------------------
! (1) Molecular Hamiltonian contribution: -i * H * v1
!
! Note that: (i)  < Psi0 | H | Psi_J > = 0 for all J (by definition
!                 of the intermediate state basis).
!
!            (ii) < Psi0 | H | Psi0 > = 0 (as we are dealing with
!                 ground state energy-shifted Hamiltonian)
!
!----------------------------------------------------------------------
! Note that here matdim is the no. intermediate states + 1
! (where the extra basis function is the ground state)
!
! Therefore the functions matxvec_treal_* will have to be passed
! matdim-1, as this is the dimension of the IS representation of the
! hamiltonian.
!
! For the same reasons, and acknowledging that the last basis function
! corresponds to the ground state, we have to pass v1(1:matdim-1)
! and v2(1:matdim-1) to the functions matxvec_treal_*.
!
! noffdiag, however, is fine as is
!----------------------------------------------------------------------

    if (hincore) then
       call matxvec_treal_incore(matdim-1,noffdiag,v1(1:matdim-1),&
            opxv1(1:matdim-1))
    else
       call matxvec_treal_ext(matdim-1,noffdiag,v1(1:matdim-1),&
            opxv1(1:matdim-1))
    endif
    
    ! Contribution to v2=dt|Psi>
    v2=v2+opxv1
    
!----------------------------------------------------------------------
! (2) Dipole-laser contribution: -i * -mu.E(t) * v1 = +i * mu.E(t) * v1
!----------------------------------------------------------------------
! Three pieces: (a) IS representation of the dipole operator, D_IJ.
!               (b) Ground state dipole matrix element, D_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, D_0J.
!----------------------------------------------------------------------
    ! External electric field
    Et=efield(time)
    
    ! (a) IS-IS block
    !
    ! Loop over components of the dipole operator
    do i=1,3
       
       ! Cycle if the  current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Calculate Dc*v1
       filename='SCRATCH/dipole_'//acomp(i)
       call opxvec_ext(matdim-1,v1(1:matdim-1),&
            opxv1(1:matdim-1),filename,nbuf_dip(i))
       
       ! Contribution to v2=dt|Psi>
       v2(1:matdim-1)=v2(1:matdim-1) &
            + ci*Et(i)*opxv1(1:matdim-1) &
            + ci*Et(i)*d00(i)*v1(1:matdim-1)

    enddo

    ! (b) Ground state-ground state element
    !
    ! Loop over components of the dipole operator
    do i=1,3
       
       ! Cycle if the  current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Contribution to v2=dt|Psi>
       v2(matdim)=v2(matdim)+ci*Et(i)*d00(i)*v1(matdim)

    enddo

    ! (c) Ground state-IS block
    !
    ! Loop over components of the dipole operator
    do i=1,3

       ! Cycle if the current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Contribution to v2=dt|Psi>
       v2(matdim)=v2(matdim)+&
            ci*Et(i)*dot_product(d0j(i,1:matdim-1),v1(1:matdim-1))
       v2(1:matdim-1)=v2(1:matdim-1)+&
            ci*Et(i)*d0j(i,1:matdim-1)*v1(matdim)

    enddo
    
!----------------------------------------------------------------------
! (3) CAP contribution: -i * -i * W * v1 = - W * v1
!----------------------------------------------------------------------
! Three pieces: (a) IS representation of the CAP operator, W_IJ.
!               (b) Ground state CAP matrix element, W_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, W_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    if (lcap) then

       ! (a) IS-IS block
       !
       ! Calculate W*v1
       filename='SCRATCH/cap'
       call opxvec_ext(matdim-1,v1(1:matdim-1),&
            opxv1(1:matdim-1),filename,nbuf_cap)
       
       ! Contribution to v2=dt|Psi>
       v2(1:matdim-1)=v2(1:matdim-1) &
            -opxv1(1:matdim-1) &
            -w00*v1(1:matdim-1)
       
       ! (b) Ground state-ground state element
       !
       ! Contribution to v2=dt|Psi>
       if (.not.lprojcap.or.statenumber.gt.0) then
          v2(matdim)=v2(matdim)-w00*v1(matdim)
       endif

       ! (c) Ground state-IS block
       !
       ! Contribution to v2=dt|Psi>
       if (.not.lprojcap.or.statenumber.gt.0) then
          v2(matdim)=v2(matdim) &
               -dot_product(w0j(1:matdim-1),v1(1:matdim-1))
          v2(1:matdim-1)=v2(1:matdim-1) &
               -w0j(1:matdim-1)*v1(matdim)

       endif

    endif
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(opxv1)
    
    return
    
  end subroutine matxvec_treal_laser

!#######################################################################
! efield: calculates the value of the applied external electric field
!         at the time t
!#######################################################################

  function efield(t)

    use constants
    use parameters
    
    implicit none

    real(d), dimension(3) :: efield
    real(d)               :: pulse,envelope
    real(d)               :: t

!----------------------------------------------------------------------
! Carrier wave value
!----------------------------------------------------------------------
    if (ipulse.eq.1) then
       ! cosine pulse
       pulse=cos(freq*(t-t0)-phase)
    else if (ipulse.eq.2) then
       ! sine pulse
       pulse=sin(freq*(t-t0)-phase)
    endif

!----------------------------------------------------------------------
! Envelope value
!----------------------------------------------------------------------
    if (ienvelope.eq.1) then
       ! Cosine squared envelope
       envelope=cos2_envelope(t)
    else if (ienvelope.eq.2) then
       ! Sine squared-ramp envelope
       envelope=sin2ramp_envelope(t)
    else if (ienvelope.eq.3) then
       ! Box-type envelope
       envelope=box_envelope(t)
    else if (ienvelope.eq.4) then
       ! Sine squared envelope
       envelope=sin2_envelope(t)
    endif
 
!----------------------------------------------------------------------
! Electric field value
!----------------------------------------------------------------------
    efield(1:3)=pulse_vec(1:3)*strength*envelope*pulse
    
    return

  end function efield

!#######################################################################
! cos2_envelope: calculates the value of a cosine squared envelope
!                function at the time t
!#######################################################################
  
  function cos2_envelope(t) result(func)

    use constants
    use parameters

    implicit none

    real(d) :: t,func,tzero,fwhm

    ! tzero
    tzero=envpar(1)

    ! fwhm
    fwhm=envpar(2)

    ! Envelope value
    if (t.gt.tzero-fwhm.and.t.lt.tzero+fwhm) then
       func=(cos(pi*(t-tzero)/(2.0d0*fwhm)))**2
    else
       func=0.0d0
    endif
    
    return
    
  end function cos2_envelope

!#######################################################################
! sin2_envelope: calculates the value of a sine squared envelope
!                function at the time t
!#######################################################################
  
  function sin2_envelope(t) result(func)

    use constants
    use parameters

    implicit none

    real(d) :: t,func,tzero,fwhm

    ! tzero
    tzero=envpar(1)

    ! fwhm
    fwhm=envpar(2)

    ! Envelope value
    if (t.gt.tzero.and.t.lt.tzero+2.0d0*fwhm) then
       func=(sin(pi*(t-tzero)/(2.0d0*fwhm)))**2
    else
       func=0.0d0
    endif
    
    return
    
  end function sin2_envelope
  
!#######################################################################
! sin2ramp_envelope: calculates the value of a sine-squared ramp
!                    envelope function at the time t.
!#######################################################################
  
  function sin2ramp_envelope(t) result(func)

    use constants
    use parameters

    implicit none

    real(d) :: t,func,tau

    ! tau
    tau=envpar(1)

    ! Envelope value
    if (t.lt.0.0d0) then
       func=0.0d0
    else if (t.ge.0.0d0.and.t.le.tau) then
       func=sin(pi*t/(2.0d0*tau))
       func=func**2
    else
       func=1.0d0
    endif

    return
    
  end function sin2ramp_envelope

!#######################################################################
! box_envelope: calculates the value of a box function at the time t.
!#######################################################################

  function box_envelope(t) result(func)

    use constants
    use parameters

    implicit none

    real(d) :: t,func,ti,tf

    ! t_i
    ti=envpar(1)

    ! t_f
    tf=envpar(2)

    ! Envelope value
    if (t.lt.ti.or.t.gt.tf) then
       func=0.0d0
    else
       func=1.0d0
    endif
    
    return
    
  end function box_envelope
    
!#######################################################################
! opxvec_treal_ext: calculates v2 = O * v1 using an externally stored
!                   operator matrix O
!#######################################################################

  subroutine opxvec_ext(matdim,v1,v2,filename,nrec)

    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
    
    implicit none

    integer                                 :: matdim,nrec,nlim,&
                                               unit,i,k,n,nthreads,&
                                               tid,npt,tmp
    integer, dimension(:), allocatable      :: indxi,indxj
    integer*8, dimension(:,:), allocatable  :: irange
    real(d), dimension(:), allocatable      :: oij,oii
    complex(d), dimension(matdim)           :: v1,v2
    complex(d), dimension(:,:), allocatable :: tmpvec
    character(len=70)                       :: filename
    
!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(irange(nthreads,2))
    allocate(tmpvec(matdim,nthreads))
    allocate(oii(matdim))
    allocate(oij(buf_size))
    allocate(indxi(buf_size))
    allocate(indxj(buf_size))
    
!-----------------------------------------------------------------------
! Open the operator matrix file
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,status='old',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements to the matrix-vector
! product
!-----------------------------------------------------------------------
    ! Read past the buffer size, which we already know
    read(unit) tmp

    ! Read the on-diagonal elements
    read(unit) oii

    ! Contribution from the on-diagonal elements
    do k=1,matdim
       v2(k)=oii(k)*v1(k)
    enddo
    
!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements to the matrix-vector
! product
!-----------------------------------------------------------------------
    ! Initialisation
    tmpvec=czero

    ! Loop over records
    do n=1,nrec

       ! Read the current record
       read(unit) oij(:),indxi(:),indxj(:),nlim
       
       ! Partitioning of the operator matrix elements: one chunk per
       ! thread
       npt=int(floor(real(nlim)/real(nthreads)))       
       do i=1,nthreads-1
          irange(i,1)=(i-1)*npt+1
          irange(i,2)=i*npt
       enddo
       irange(nthreads,1)=(nthreads-1)*npt+1
       irange(nthreads,2)=nlim

       ! Contribution of the current record to the matrix-vector
       ! product
       !
       !$omp parallel do &
       !$omp& private(i,k,tid) &
       !$omp& shared(tmpvec,v1,oij,indxi,indxj)
       do i=1,nthreads
          tid=1+omp_get_thread_num()
          do k=irange(tid,1),irange(tid,2)
             tmpvec(indxi(k),tid)=tmpvec(indxi(k),tid)&
                  +oij(k)*v1(indxj(k))
             tmpvec(indxj(k),tid)=tmpvec(indxj(k),tid)&
                  +oij(k)*v1(indxi(k))
          enddo
       enddo
       !$omp end parallel do

    enddo

    do i=1,nthreads
       v2(:)=v2(:)+tmpvec(:,i)
    enddo
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(irange)
    deallocate(tmpvec)
    deallocate(oii)
    deallocate(oij)
    deallocate(indxi)
    deallocate(indxj)
    
!-----------------------------------------------------------------------
! Close the operator matrix file
!-----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine opxvec_ext

!#######################################################################

  subroutine matxvec_treal_laser_adc1(time,matdim,noffdiag,v1,v2)

    use constants
    use parameters
    
    implicit none
    
    integer                       :: matdim,i,j
    integer*8                     :: noffdiag
    real(d)                       :: time
    real(d), dimension(3)         :: Et
    complex(d), dimension(matdim) :: v1,v2

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    v2=czero

!----------------------------------------------------------------------
! (1) Molecular Hamiltonian contribution: -i * H * v1
!
! Note that: (i)  < Psi0 | H | Psi_J > = 0 for all J (by definition
!                 of the intermediate state basis).
!
!            (ii) < Psi0 | H | Psi0 > = 0 (as we are dealing with
!                 ground state energy-shifted Hamiltonian)
!
!----------------------------------------------------------------------
    v2(1:matdim-1)=v2(matdim-1)-ci*matmul(h1,v1(1:matdim-1))

!----------------------------------------------------------------------
! (2) Dipole-laser contribution: -i * -mu.E(t) * v1 = +i * mu.E(t) * v1
!----------------------------------------------------------------------
! Three pieces: (a) IS representation of the dipole operator, D_IJ.
!               (b) Ground state dipole matrix element, D_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, D_0J.
!----------------------------------------------------------------------
    ! External electric field
    Et=efield(time)

    ! (a) IS-IS block
    !
    ! Loop over components of the dipole operator
    do i=1,3
       
       ! Cycle if the current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Contribution to v2=dt|Psi>
       v2(1:matdim-1)=v2(1:matdim-1)&
            + ci*Et(i)*matmul(dij(i,:,:),v1(1:matdim-1)) &
            + ci*Et(i)*d00(i)*v1(1:matdim-1)

    enddo

    ! (b) Ground state-ground state element
    !
    ! Loop over components of the dipole operator
    do i=1,3
       
       ! Cycle if the  current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Contribution to v2=dt|Psi>
       v2(matdim)=v2(matdim)+ci*Et(i)*d00(i)*v1(matdim)

    enddo

    ! (c) Ground state-IS block
    !
    ! Loop over components of the dipole operator
    do i=1,3

       ! Cycle if the  current component does not contribute
       if (pulse_vec(i).eq.0.0d0) cycle

       ! Contribution to v2=dt|Psi>
       v2(matdim)=v2(matdim)+&
            ci*Et(i)*dot_product(d0j(i,1:matdim-1),v1(1:matdim-1))
       v2(1:matdim-1)=v2(1:matdim-1)+&
            ci*Et(i)*d0j(i,1:matdim-1)*v1(matdim)
       
    enddo

!----------------------------------------------------------------------
! (3) CAP contribution: -i * -i * W * v1 = - W * v1
!----------------------------------------------------------------------
! Three pieces: (a) IS representation of the CAP operator, W_IJ.
!               (b) Ground state CAP matrix element, W_00.
!               (c) Off-diagonal elements between the ground state
!                   and the intermediate states, W_0J.
!
! Note that (b) and (c) do not contribute if the CAP is projected
! onto the space orthogonal to the ground state
!----------------------------------------------------------------------
    if (lcap) then

       ! (a) IS-IS block
       !
       v2(1:matdim-1)=v2(1:matdim-1) &
            -matmul(wij,v1(1:matdim-1)) &
            -w00*v1(1:matdim-1)
       
       ! (b) Ground state-ground state element
       !
       if (.not.lprojcap.or.statenumber.gt.0) then
          v2(matdim)=v2(matdim)-w00*v1(matdim)
       endif

       ! (c) Ground state-IS block
       !
       if (.not.lprojcap.or.statenumber.gt.0) then
          v2(matdim)=v2(matdim) &
               -dot_product(w0j(1:matdim-1),v1(1:matdim-1))
          v2(1:matdim-1)=v2(1:matdim-1) &
               -w0j(1:matdim-1)*v1(matdim)
       endif

    endif
    
    return
    
  end subroutine matxvec_treal_laser_adc1
    
!#######################################################################
! chebyshev_recursion: Calculates the Chebyshev vectors using the
!                      two-term Chebyshev recursion relation.
!
!                      Used in Chebyshev propagation and the calculation
!                      of the Chebyshev order-domain autocorrelation
!                      function.
!
!                      It is debatable whether this is the most suitable
!                      module in which to have this routine, but putting
!                      it here allows us to make use of the
!                      load_hamiltonian, etc routines.
!#######################################################################  

  subroutine chebyshev_recursion(k,matdim,noffdiag,bounds,qk,qkm1,qkm2)

    use constants
    
    implicit none

    integer, intent(in)               :: k,matdim
    integer*8, intent(in)             :: noffdiag
    real(d), dimension(2), intent(in) :: bounds
    real(d), dimension(matdim)        :: qk,qkm1,qkm2

!----------------------------------------------------------------------
! 2 * H_norm * q_k-1 (k>2)
! or  H_norm * q_k-1 (k=1)
!----------------------------------------------------------------------
    ! H_norm * q_k-1
    call matxvec_chebyshev(matdim,noffdiag,bounds,qkm1,qk)

    ! 2* H_norm * q_k-1
    if (k.gt.1) qk=2.0d0*qk
       
!----------------------------------------------------------------------
! q_k = 2 * H_norm * q_k-1 - q_k-2 (k>2)
!----------------------------------------------------------------------
    if (k.gt.1) qk=qk-qkm2

    return
    
  end subroutine chebyshev_recursion

!#######################################################################
! matxvec_chebyshev: Calculates v2 = H_norm * v1, where H_norm is the
!                    Hamiltonian matrix scaled to have eigenvalues in
!                    the closed interval [-1,1]
!#######################################################################
  
  subroutine matxvec_chebyshev(matdim,noffdiag,bounds,v1,v2)
    
    use constants
    use iomod
    
    implicit none

    integer, intent(in)               :: matdim
    integer*8, intent(in)             :: noffdiag
    real(d), dimension(2), intent(in) :: bounds
    real(d), dimension(matdim)        :: v1,v2

    if (hincore) then
       call matxvec_chebyshev_incore(matdim,noffdiag,bounds,v1,v2)
    else
       errmsg='WRITE THE OUT-OF-CORE MATXVEC_CHEBYSHEV CODE'
       call error_control
    endif
    
    return
    
  end subroutine matxvec_chebyshev

!#######################################################################
! matxvec_chebyshev_incore: In-core calculation of v2 = +H_norm * v1,
!                           where H_norm is the Hamiltonian matrix
!                           scaled to have eigenvalues in the closed
!                           interval [-1,1]
!#######################################################################

  subroutine matxvec_chebyshev_incore(matdim,noffdiag,bounds,v1,v2)
    
    use constants
    use parameters
    use iomod
    use channels
    use omp_lib
    
    implicit none

    integer, intent(in)                    :: matdim
    integer*8, intent(in)                  :: noffdiag
    integer                                :: nthreads,i,k,tid
    integer*8, allocatable                 :: irange(:,:)
    integer*8                              :: npt
    real(d), dimension(2), intent(in)      :: bounds
    real(d), dimension(matdim), intent(in) :: v1
    real(d), dimension(matdim)             :: v2
    real(d), allocatable                   :: tmpvec(:,:)
    real(d)                                :: DeltaE,Emin
    
!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(irange(nthreads,2))
    allocate(tmpvec(matdim,nthreads))

!-----------------------------------------------------------------------
! Partitioning of the off-diagonal elements: one chunk per thread
!-----------------------------------------------------------------------
    npt=int(floor(real(noffdiag)/real(nthreads)))
    
    do i=1,nthreads-1
       irange(i,1)=(i-1)*npt+1
       irange(i,2)=i*npt
    enddo

    irange(nthreads,1)=(nthreads-1)*npt+1
    irange(nthreads,2)=noffdiag

!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements to H*v1
!-----------------------------------------------------------------------
    do k=1,matdim
       v2(k)=hon(k)*v1(k)
    enddo

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements to H*v1
!-----------------------------------------------------------------------
    tmpvec=0.0d0

    !$omp parallel do &
    !$omp& private(i,k,tid) &
    !$omp& shared(tmpvec,v1)
    do i=1,nthreads
       
       tid=1+omp_get_thread_num()

       do k=irange(tid,1),irange(tid,2)
          tmpvec(iindx(k),tid)=tmpvec(iindx(k),tid)&
               +hoff(k)*v1(jindx(k))
          tmpvec(jindx(k),tid)=tmpvec(jindx(k),tid)&
               +hoff(k)*v1(iindx(k))
       enddo
          
    enddo
    !$omp end parallel do

    do i=1,nthreads
       v2=v2+tmpvec(:,i)
    enddo

!-----------------------------------------------------------------------
! DeltaE and Emin
!-----------------------------------------------------------------------
    DeltaE=bounds(2)-bounds(1)
    Emin=bounds(1)
    
!-----------------------------------------------------------------------
! (2/DeltaE)*H*v1
!-----------------------------------------------------------------------
    v2=(2.0d0/DeltaE)*v2

!-----------------------------------------------------------------------
! (2/DeltaE)*H*v1 - 2*[(0.5*DeltaE + E_min)/DeltaE]*v1
!-----------------------------------------------------------------------
    v2=v2-2.0d0*(0.5d0*DeltaE+Emin)*(1.0d0/DeltaE)*v1
    
!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(irange)
    deallocate(tmpvec)
    
    return

  end subroutine matxvec_chebyshev_incore
  
!#######################################################################
  
end module tdsemod
