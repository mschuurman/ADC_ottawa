!#######################################################################
! mp2: routines for the calculation of various MP2 quantities: energies
!      density matrices and the D2 diagnostic
!#######################################################################

module mp2

contains

!#######################################################################

  subroutine mp2_master(e0)

    use channels
    use constants
    use parameters
    
    implicit none
        
    integer           :: i
    real(dp)          :: e0,d2
    character(len=60) :: atmp
    
!-----------------------------------------------------------------------
! Calculation of the MP2 correlation energy
!-----------------------------------------------------------------------
    call mp2_ener

!-----------------------------------------------------------------------
! Calculation of the D2 diagnostic
!-----------------------------------------------------------------------
    call mp2_d2(d2)

!-----------------------------------------------------------------------
! Output results
!-----------------------------------------------------------------------
    e0=Ehf+E_MP2
 
    write(ilog,'(/,2x,90a)') ('*',i=1,70)
    
    write(ilog,'(2x,1a)') '*'
    atmp='* HF energy:'
    write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),Ehf
    
    write(ilog,'(2x,1a)') '*'
    atmp='* MP2 correlation energy:'
    write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),E_MP2
    
    write(ilog,'(2x,1a)') '*'
    atmp='* MP2 energy:'
    write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),e0
    
    write(ilog,'(2x,1a)') '*'
    atmp='* D2 diagnostic:'
    write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),d2
    
    write(ilog,'(2x,1a)') '*'
    write(ilog,'(2x,90a,/)') ('*',i=1,70)

    return

  end subroutine mp2_master

!#######################################################################
! mp2_ener: Calculates the MP2 correlation energy
!#######################################################################

  subroutine mp2_ener

    use channels
    use constants
    use parameters
    use vpqrsmod
    use omp_lib

    implicit none

    integer                             :: r,s,u,v,nsym1,nsym2,i,j,a,b
    integer                             :: occ1,nthreads,tid
    real(dp)                            :: DA,eijc,term,etotal
    real(dp), dimension(:), allocatable :: sum1thread

!-----------------------------------------------------------------------
! Set the index of the first non-frozen occupied orbital
!-----------------------------------------------------------------------
    if (lifrzcore) then
       occ1=ncore+1
    else
       occ1=1
    endif

!-----------------------------------------------------------------------    
! Determine the number of threads
!-----------------------------------------------------------------------    
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    allocate(sum1thread(nthreads))
    sum1thread=0.0d0

!-----------------------------------------------------------------------    
! Calculate the MP2 energy correction
!-----------------------------------------------------------------------    
    E_MP2 = 0.0d0
    etotal = 0.0d0

    !$omp parallel do &
    !$omp& private(r,a,s,b,u,i,v,j,nsym1,nsym2,eijc,term,tid) &
    !$omp& shared(sum1thread,e,orbsym,roccnum,MT)
    do r=nOcc+1,nBas
       a=roccnum(r)  !r
           
       do s=nOcc+1,nBas
          b=roccnum(s) !s
          
          do u=occ1,nocc
             i=roccnum(u) !u
                 
             do v=occ1,nocc
                j=roccnum(v) !v

                tid=1+omp_get_thread_num()
                
                nsym1=MT(orbSym(i),orbSym(a))
                nsym2=MT(orbSym(j),orbSym(b))
                    
                if  (MT(nsym1,nsym2) .eq. 1)  then

                   eijc=e(i)+e(j)-e(a)-e(b)
                       
                   term= vpqrs(i,a,j,b)*(2._dp*vpqrs(i,a,j,b)-vpqrs(i,b,j,a))
                       
                   term=term/eijc

                   sum1thread(tid)=sum1thread(tid)+term
                       
                endif
                    
             enddo

          enddo

       enddo

    enddo
    !$omp end parallel do

    E_MP2=sum(sum1thread)
    
    deallocate(sum1thread)
    
    return
    
  end subroutine mp2_ener

!#######################################################################
! mp2_d2: Calculates the D2(MP1) MP2 diagnostic of Nielsen and
!         Janssen. See Chem. Phys. Lett., 310, 568 (1999).
!#######################################################################

  subroutine mp2_d2(d2)

    use channels
    use constants
    use parameters
    use iomod
    use omp_lib
      
    implicit none

    integer                               :: i,j,k,a,b,c,albl,blbl,&
                                             workdim,error
    real(dp)                              :: d2
    real(dp), dimension(:,:), allocatable :: toto,tvtv
    real(dp), dimension(:), allocatable   :: eigval,workarr

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(toto(nocc,nocc))
    allocate(tvtv(nvirt,nvirt))

!-----------------------------------------------------------------------
! Calculation of T^o * T^o
!-----------------------------------------------------------------------
    toto=0.0d0

    !$omp parallel do &
    !$omp& private(i,j,k,a,b) &
    !$omp& shared(toto)
    do i=1,nocc
       do j=1,nocc
          do k=1,nocc
             do a=nocc+1,nbas
                do b=nocc+1,nbas
                   toto(i,j)=toto(i,j)+t2amp(i,k,a,b)*t2amp(j,k,a,b)
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
      
!-----------------------------------------------------------------------
! Calculation of T^v * T^v
!-----------------------------------------------------------------------
    tvtv=0.0d0
      
    !$omp parallel do &
    !$omp& private(a,albl,b,blbl,i,j,c) &
    !$omp& shared(tvtv)
    do a=nocc+1,nbas
       albl=a-nocc
       do b=nocc+1,nbas
          blbl=b-nocc            
          do i=1,nocc
             do j=1,nocc
                do c=nocc+1,nbas
                   tvtv(albl,blbl)=tvtv(albl,blbl)+&
                        t2amp(i,j,a,c)*t2amp(i,j,b,c)
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
      
!-----------------------------------------------------------------------
! Eigenvalues of T^o * T^o
!-----------------------------------------------------------------------
    workdim=3*nocc
    allocate(workarr(workdim))
    allocate(eigval(nocc))
    
    call dsyev('N','U',nocc,toto,nocc,eigval,workarr,workdim,error)
    
    if (error.ne.0) then
       errmsg='Diagonalisation of T^o * T^o failed in subroutine &
            mp2_d2'
       call error_control
    endif
    
    d2=sqrt(eigval(nocc))
    
    deallocate(workarr)
    deallocate(eigval)

!-----------------------------------------------------------------------
! Eigenvalues of T^v * T^v
!-----------------------------------------------------------------------
    workdim=3*nvirt
    allocate(workarr(workdim))
    allocate(eigval(nvirt))

    call dsyev('N','U',nvirt,tvtv,nvirt,eigval,workarr,workdim,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of T^v * T^v failed in subroutine &
            mp2_d2'
       call error_control
    endif

    if (sqrt(eigval(nvirt)).gt.d2) d2=sqrt(eigval(nvirt))
    
    deallocate(workarr)
    deallocate(eigval)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(toto)
    deallocate(tvtv)

    return

  end subroutine mp2_d2

!#######################################################################

  function t2amp(i,j,a,b)
    
    use constants
    use parameters
    use vpqrsmod
    
    implicit none
    
    integer, intent(in)  :: i,j,a,b
    real(dp)             :: t2amp
    
    t2amp=vpqrs(i,a,j,b)/(e(i)+e(j)-e(a)-e(b))
    
    return
    
  end function t2amp

!#######################################################################
! rho_mp2: calculates the MP2 ground state density matrix
!#######################################################################

  subroutine rho_mp2(rho)

    use channels
    use constants
    use parameters
    use timingmod
    use iomod

    implicit none

    integer                        :: i,j,a,b,k,occ1
    real(dp), dimension(nbas,nbas) :: rho
    real(dp)                       :: tw1,tw2,tc1,tc2

    integer                        :: error,workdim
    real(dp), allocatable          :: natvec(:,:)
    real(dp), allocatable          :: natocc(:)
    real(dp), allocatable          :: work(:)
    
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)

!-----------------------------------------------------------------------
! Ouput what we are doing
!-----------------------------------------------------------------------
    write(ilog,'(2x,a)') 'Calculating the MP2 density matrix'

!-----------------------------------------------------------------------
! Set the index of the first non-frozen occupied orbital
!-----------------------------------------------------------------------
    if (lifrzcore) then
       occ1=ncore+1
    else
       occ1=1
    endif

!-----------------------------------------------------------------------
! Ground state density matrix.
! Note that the 1st-order correction is zero.
!-----------------------------------------------------------------------
    rho=0.0d0

    ! Occupied-occupied block: 0th-order contribution
    do i=1,nocc
       rho(i,i)=2.0d0
    enddo

    ! Occupied-occupied block: 2nd-order contribution
    do i=occ1,nocc
       do j=occ1,nocc
          rho(i,j)=rho(i,j)+rhogs2_oo(i,j)
          rho(j,i)=rho(i,j)
       enddo
    enddo
    
    ! Unoccupied-unoccupied block: 2nd-order contribution
    do a=nocc+1,nbas
       do b=a,nbas
          rho(a,b)=rhogs2_uu(a,b)
          rho(b,a)=rho(a,b)
       enddo
    enddo
    
    ! Occupied-unoccupied block: 2nd-order contribution
    do i=occ1,nocc
       do a=nocc+1,nbas
          rho(i,a)=rhogs2_ou(i,a)
          rho(a,i)=rho(i,a)
       enddo
    enddo

!-----------------------------------------------------------------------
! Check: natural orbital occupation numbers
!-----------------------------------------------------------------------
    !workdim=3*nbas
    !allocate(work(workdim))
    !allocate(natvec(nbas,nbas))
    !allocate(natocc(nbas))
    !
    !natvec=rho
    !
    !call dsyev('N','U',nbas,natvec,nbas,natocc,work,workdim,error)
    !
    !if (error.ne.0) then
    !   errmsg='Error in the diagonalisation of the MP2 density matrix'
    !endif
    !
    !print*,
    !do i=1,nbas
    !   write(6,'(i2,2x,F7.4)') i,natocc(nbas-i+1)
    !enddo
    !print*,
    !print*,sum(natocc(1:nbas))
    !print*,sum(natocc(nbas-nocc+1:nbas))
    !print*,sum(natocc(1:nbas-nocc))
    !print*,
    !STOP
    !
    !deallocate(work)
    !deallocate(natvec)
    !deallocate(natocc)

!-----------------------------------------------------------------------
! Finish timing and output the wall time taken
!-----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(2x,a,2x,F7.2,1x,a1,/)') 'Time taken:',tw2-tw1,'s'
    
    return
    
  end subroutine rho_mp2
  
!#######################################################################
! rhogs2_uu: returns an element of the unoccupied-unoccupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################
   
  function rhogs2_uu(a,b) result(fret)

    use constants
    use parameters
    use vpqrsmod
    
    implicit none

    integer  :: a,b,i1,j1,c1,i,j,c,occ1
    real(dp) :: fret,ftmp,delta_ijac,delta_ijbc
      
!-----------------------------------------------------------------------
! Set the index of the first non-frozen occupied orbital
!-----------------------------------------------------------------------
    if (lifrzcore) then
       occ1=ncore+1
    else
       occ1=1
    endif

!-----------------------------------------------------------------------
! Calculate the matrix element
!-----------------------------------------------------------------------
    fret=0.0d0

    do i1=occ1,nocc
       i=roccnum(i1)
       
       do j1=occ1,nocc
          j=roccnum(j1)
       
          do c1=nocc+1,nbas
             c=roccnum(c1)
             
             delta_ijac=1.0d0/(e(a)+e(c)-e(i)-e(j))
             delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))
                          
             ftmp=vpqrs(a,i,c,j)*(2.0d0*vpqrs(i,b,j,c)-vpqrs(i,c,j,b))&
                  +vpqrs(c,i,a,j)*(2.0d0*vpqrs(i,c,j,b)-vpqrs(i,b,j,c))
             
             fret=fret+delta_ijac*delta_ijbc*ftmp
             
          enddo
       enddo
    enddo

    return

  end function rhogs2_uu

!#######################################################################
! rhogs2_oo: returns an element of the occupied-occupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################

  function rhogs2_oo(i,j) result(fret)

    use constants
    use parameters
    use vpqrsmod
    
    implicit none

    integer  :: i,j,k1,a1,b1,k,a,b,occ1
    real(dp) :: fret,ftmp,delta_ikab,delta_jkab

!-----------------------------------------------------------------------
! Set the index of the first non-frozen occupied orbital
!-----------------------------------------------------------------------
    if (lifrzcore) then
       occ1=ncore+1
    else
       occ1=1
    endif
    
!-----------------------------------------------------------------------
! Calculate the matrix element
!-----------------------------------------------------------------------
    fret=0.0d0

    do k1=occ1,nocc
      k=roccnum(k1)
       do a1=nocc+1,nbas
          a=roccnum(a1)
          do b1=nocc+1,nbas
             b=roccnum(b1)
             
             delta_ikab=1.0d0/(e(a)+e(b)-e(i)-e(k))
             delta_jkab=1.0d0/(e(a)+e(b)-e(j)-e(k))
             
             ftmp=vpqrs(a,i,b,k)*(2.0d0*vpqrs(j,a,k,b)-vpqrs(j,b,k,a))&
                  +vpqrs(b,i,a,k)*(2.0d0*vpqrs(j,b,k,a)-vpqrs(j,a,k,b))
             
             fret=fret-delta_ikab*delta_jkab*ftmp
             
          enddo
       enddo
    enddo

    return

  end function rhogs2_oo

!#######################################################################
! rhogs2_ou: returns an element of the occupied-unoccupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################

  function rhogs2_ou(i,a) result(fret)

    use constants
    use parameters
    use vpqrsmod
    
    implicit none

    integer  :: i,a,j1,k1,b1,c1,j,k,b,c,occ1
    real(dp) :: fret,ftmp,term1,term2,delta_ijbc,delta_jkab

!-----------------------------------------------------------------------
! Set the index of the first non-frozen occupied orbital
!-----------------------------------------------------------------------
    if (lifrzcore) then
       occ1=ncore+1
    else
       occ1=1
    endif

!-----------------------------------------------------------------------
! Calculate the matrix element
!-----------------------------------------------------------------------
    term1=0.0d0
    do j1=occ1,nocc
       j=roccnum(j1)
       do b1=nocc+1,nbas
          b=roccnum(b1)
          do c1=nocc+1,nbas
             c=roccnum(c1)
             
             delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))
      
             ftmp=vpqrs(i,b,j,c)*(vpqrs(j,b,a,c)-2.0d0*vpqrs(j,c,a,b))&
                  +vpqrs(i,c,j,b)*(vpqrs(j,c,a,b)-2.0d0*vpqrs(j,b,a,c))

             term1=term1+delta_ijbc*ftmp

          enddo
       enddo
    enddo

    term2=0.0d0
    do j1=occ1,nocc
       j=roccnum(j1)
       do k1=occ1,nocc
          k=roccnum(k1)
          do b1=nocc+1,nbas
             b=roccnum(b1)
             
             delta_jkab=1.0d0/(e(a)+e(b)-e(j)-e(k))
             
             ftmp=vpqrs(j,i,k,b)*(2.0d0*vpqrs(j,a,k,b)-vpqrs(j,b,k,a))&
                  +vpqrs(j,b,k,i)*(2.0d0*vpqrs(j,b,k,a)-vpqrs(j,a,k,b))
             
             term2=term2+delta_jkab*ftmp

          enddo
       enddo
    enddo

    fret=-(term1+term2)/(e(a)-e(i))
    
    return

  end function rhogs2_ou
    
!#######################################################################

  end module mp2
