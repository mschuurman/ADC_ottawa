!######################################################################
! Routines to compute the unitary transformation matrix that
! maximally approximately diagonalises a set of matrices
!######################################################################

module simdiag

contains

!######################################################################
! simdiag_jacobi: Simultaneous diagonalisation of a set of matrices
!                 using a Jacobi rotation scheme.
!
!                 For a set of real symmetric matrices A, returns an
!                 orthogonal transformation matrix Q.
!
!                 Uses the algorithms reported in
!
!                 SIAM J, Matrix Anal. Appl. 14, 927 (1993)
!
!                 and
!
!                 SIAM J, Matrix Anal. Appl. 17, 161 (1996)
!######################################################################
    
  subroutine simdiag_jacobi(A1,Q,dim,nmat)

    use constants
    use iomod
    use channels
    use timingmod
    
    implicit none

    integer                          :: dim,nmat
    integer                          :: i,j,k,l,m,n,noff,error
    integer, parameter               :: maxit=250
    real(d), dimension(dim,dim,nmat) :: A1,A
    real(d), dimension(dim,dim)      :: Q,tmp
    real(d), dimension(dim)          :: lambda
    real(d), dimension(3*dim)        :: work
    real(d)                          :: off,c,s,cs,sc
    real(d)                          :: last
    real(d), parameter               :: thrsh=1e-12_d
    real(d)                          :: tw1,tw2,tc1,tc2
    logical                          :: converged

!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(2x,a)') 'Approximate simultaneous diagonalisation &
         of X, Y, and Z'
    write(ilog,'(2x,a,/)') 'Alorithm: Jacobi'

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Initialisation: set the orthogonal transformation matrix Q as the
! matrix that diagonalises one of the sum matrices in the set A
!----------------------------------------------------------------------
    ! Number of off-diagonal elements in the lower triangle
    noff=dim*(dim-1)/2    

    ! Make a copy of the set of matrices to diagonalise
    A=A1

    ! Diagonalise sum_k A_k/nmat
    Q=0.0d0
    do k=1,nmat
       Q=Q+A(:,:,k)/nmat
    enddo
    call dsyev('V','U',dim,Q,dim,lambda,work,3*dim,error)    
    if (error.ne.0) then
       errmsg='Initialisation failure in subroutine simdiag_jacobi'
       call error_control
    endif
    
    ! Similarity transform all matrices in the set
    do k=1,nmat
       tmp=A(:,:,k)
       A(:,:,k)=matmul(transpose(Q),matmul(tmp,Q))
    enddo

    ! Compute and save the initial objective function value
    off=0.0d0
    do k=1,nmat
       do i=1,dim-1
          do j=i+1,dim
             off=off+A(j,i,k)**2/(nmat*noff)
          enddo
       enddo
    enddo
    last=off

    ! Output our progress
    write(ilog,'(2x,a10,i3,a2,ES15.7,2x,ES15.7)') &
         'Iteration ',n,': ',off,abs(last-off)
    
!----------------------------------------------------------------------
! Perform the simultaneous diagonalisation
!----------------------------------------------------------------------
    converged=.false.

    do n=1,maxit

       ! Loop over unique pairs of off-diagonal elements
       do i=1,dim-1
          do j=i+1,dim

             ! Compute cosine and sine of the rotation angle
             call get_cs(i,j,A,c,s,dim,nmat)
             
             ! Update the matrices in the set A
             do l=1,nmat

                ! Recalculate rows i and j
                do k=1,dim
                   cs=c*A(i,k,l)+s*A(j,k,l)
                   sc=-s*A(i,k,l)+c*A(j,k,l)
                   A(i,k,l)=cs
                   A(j,k,l)=sc
                enddo
                ! New matrix A_{n+1} from A_{n}
                do k=1,dim
                   cs=c*A(k,i,l)+s*A(k,j,l)
                   sc=-s*A(k,i,l)+c*A(k,j,l)
                   A(k,i,l)=cs
                   A(k,j,l)=sc
                enddo
             enddo

             ! Update the matrix of eigenvectors
             tmp=Q
             Q(:,i)=c*tmp(:,i)+s*tmp(:,j)
             Q(:,j)=c*tmp(:,j)-s*tmp(:,i)
             
          enddo
       enddo
       
       ! Compute the objective function value
       off=0.0d0
       do k=1,nmat
          do i=1,dim-1
             do j=i+1,dim
                off=off+A(j,i,k)**2/(nmat*noff)
             enddo
          enddo
       enddo
       
       ! Output our progress
       write(ilog,'(2x,a10,i3,a2,ES15.7,2x,ES15.7)') &
            'Iteration ',n,': ',off,abs(last-off)
       
       ! Exit if we have reached convergence
       if (abs(last-off).lt.thrsh) then
          converged=.true.
          exit
       else
          last=off
       endif
          
    enddo

!----------------------------------------------------------------------
! Exit if convergence was not reached
!----------------------------------------------------------------------
    if (.not.converged) then
       errmsg='Convergence not reaced in simdiag_jacobi'
       call error_control
    endif

!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F11.4,1x,a)') 'Time taken:',tw2-tw1," s"

!-----------------------------------------------------------------------    
! Check
!-----------------------------------------------------------------------    
    !do k=1,nmat
    !   A(:,:,k)=matmul(transpose(Q),matmul(A1(:,:,k),Q))
    !enddo
    !off=0.0d0
    !do k=1,nmat
    !   do i=1,dim-1
    !      do j=i+1,dim
    !         off=off+A(i,j,k)**2/(nmat*noff)
    !      enddo
    !   enddo
    !enddo
    !print*,
    !print*,'off: ',off
    !print*,
    !stop
    
    return
      
  end subroutine simdiag_jacobi

!######################################################################

  subroutine get_cs(i,j,A,c,s,dim,nmat)

    use constants
    use iomod

    integer                          :: i,j,dim,nmat
    integer                          :: k,p,q,error
    real(d)                          :: c,s
    real(d)                          :: x,y,r
    real(d), dimension(dim,dim,nmat) :: A
    real(d), dimension(2,nmat)       :: h
    real(d), dimension(2,2)          :: G,eigvec
    real(d), dimension(2)            :: eigval
    real(d), dimension(6)            :: work
    real(d)                          :: tiny

    real(d) :: theta,val
    real(d) :: t_on,t_off
  
!----------------------------------------------------------------------
! Cutoff for taking double precision values to be zero
!----------------------------------------------------------------------
    tiny=epsilon(tiny)

!----------------------------------------------------------------------
! TEST: Jacobi angle for a single matrix
!----------------------------------------------------------------------
    !! TEST: Jacobi angle for a single matrix
    !if (A(j,j,1)-A(i,i,1).ne.0.0d0) then
    !   val=2.0d0*A(i,j,1)/(A(j,j,1)-A(i,i,1))
    !   theta=atan(val)/2.0d0
    !else
    !   theta=pi/4.0d0
    !endif
    !c=cos(theta)
    !s=-sin(theta)
    !return
    !! TEST: Jacobi angle for a single matrix

!----------------------------------------------------------------------
! Calculate the h-vectors
!----------------------------------------------------------------------
    do k=1,nmat
       h(1,k)=A(i,i,k)-A(j,j,k)
       h(2,k)=2.0d0*A(i,j,k)
    enddo

!----------------------------------------------------------------------
! Calculate the G-matrix
!----------------------------------------------------------------------
    G=0.0d0
    do k=1,nmat
       do p=1,2
          do q=1,2
             G(p,q)=G(p,q)+h(p,k)*h(q,k)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Diagonalise the G-matrix
!----------------------------------------------------------------------
    eigvec=G

    call dsyev('V','U',2,eigvec,2,eigval,work,6,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the G-matrix failed in subroutine &
            get_cs'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate the cosine and sine of the rotation angle
! Note that: (1) We only consider inner rotations, which can be
!                achieved by ensuring that x>=0
!            (2) The constraint c**2 + s**2 = 1 has to hold
!----------------------------------------------------------------------    
!    print*,

    if (abs(eigvec(1,2)).ge.tiny) then
       !print*,"CASE I"
       x=abs(eigvec(1,2))
       y=eigvec(2,2)

       ! CHECK
       if (y.lt.1.0d0/sqrt(2.0d0)) y=-y
       ! CHECK

       r=sqrt(x**2+y**2)
       c=sqrt((x+r)/(2.0d0*r))
       s=y/(sqrt(2.0d0*r*(x+r)))
    else if (abs(eigvec(1,1)).ge.tiny) then
       !print*,"CASE II"
       x=abs(eigvec(1,1))
       y=eigvec(2,1)

       ! CHECK
       if (y.lt.1.0d0/sqrt(2.0d0)) y=-y
       ! CHECK

       r=sqrt(x**2+y**2)
       c=sqrt((x+r)/(2.0d0*r))
       s=y/(sqrt(2.0d0*r*(x+r)))
    else
       !print*,"CASE III"
       c=cos(pi/4.0d0)
       s=sin(pi/4.0d0)
    endif

    !! TEST: Single matrix case
    !print*,'ACTUAL: ',c,s
    !if (A(j,j,1)-A(i,i,1).ne.0.0d0) then
    !   val=2.0d0*A(i,j,1)/(A(j,j,1)-A(i,i,1))
    !   theta=atan(val)/2.0d0
    !else
    !   theta=pi/4.0d0
    !endif
    !!c=cos(theta)
    !!s=-sin(theta)
    !print*,'CORRECT:',cos(theta),-sin(theta)
    !! TEST: Single matrix case

    return
    
  end subroutine get_cs

!######################################################################

  function arctan2(x,y)

    use constants
    
    implicit none

    real(d) :: x,y,arctan2

    if (x.gt.0.0d0) then
       arctan2=atan(x/y)
    else if (x.lt.0.0d0.and.y.ge.0.0d0) then
       arctan2=atan(x/y)+pi
    else if (x.lt.0.0d0.and.y.lt.0.0d0) then
       arctan2=atan(x/y)-pi
    else if (x.eq.0.0d0.and.y.gt.0) then
       arctan2=0.5d0*pi
    else if (x.eq.0.0d0.and.y.lt.0) then
       arctan2=-0.5d0*pi
    else
       ! Undefined: x=0 and y=0, but compute this just so that NaN
       ! is returned
       arctan2=atan(x/y)
    endif
       
    return
    
  end function arctan2
    
!######################################################################
    
end module simdiag
